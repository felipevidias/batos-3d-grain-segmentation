"""BA-TOS v3 (scale-free) single-crop driver; also runs the v2 fixflat preset through the same code path.

Stage order and inputs are those of the fixflat harness (batos_v2.sensitivity.CropRun.execute):
  scale -> openrec -> flatten (batos_v2.flatten_fixed --variant A) -> legacy ToS seed -> background -> legacy markers
  -> true_dual K -> e10 spatialize -> e11 ToS-MIN markers -> EDT-Otsu control -> fixed-mask watershed (ToS route)
  -> E23c route -> fixed-mask watershed (selective route).
Frozen scripts are called unchanged (batos_v2.frozen.SCRIPTS); the three scripts with hard-coded constants are called
as the generated copies in batos_v3/patched (defaults = frozen values).

Presets
  fixflat : e13 scale (params JSON from the frozen e13 run), batos_v2.sensitivity.BASE_FLAGS, flatten --variant A.
            Must reproduce batos-paper-v2-fixflat bit-for-bit (checked on dev crops).
  v3      : r_granulo scale (batos_v3.scale), s = r_granulo / R_REF, flags from batos_v3.constants.stage_flags(s),
            stage-03 dark threshold = Otsu of the selector on specimen voxels (--dark-threshold), and
            "no MIN candidate inside internal_dark" recorded as an empty result instead of a crash.
GT is never read here.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import time
from pathlib import Path

import numpy as np
import tifffile

from batos_v2.extinction_run import run_tracked
from batos_v2.frozen import SCRIPTS
from batos_v2.run_crop import git_state
from batos_v2.sensitivity import BASE_FLAGS
from batos_v3 import constants
from batos_v3.scale import otsu_u8

PY = sys.executable
PATCHED = Path(__file__).resolve().parent / "patched"
SCRIPT = {
    "openrec": SCRIPTS["openrec_v2"], "openrec_lowmem": SCRIPTS["openrec_lowmem"],
    "tos_seed": SCRIPTS["batos_tos_seed"], "background": SCRIPTS["batos_background"],
    "markers_legacy": SCRIPTS["batos_markers"], "true_dual": SCRIPTS["true_dual"], "spatialize": SCRIPTS["spatialize"],
    "markers_tos": SCRIPTS["e11_markers"], "edt_otsu": PATCHED / "e11_build_control_markers_v3.py",
    "ws_tos": PATCHED / "e11_run_fixed_mask_watershed_v3.py", "ws_sel": PATCHED / "e11_run_fixed_mask_watershed_v3.py",
    "route": PATCHED / "e23c_build_selective_marker_route_v3.py",
}
LOWMEM_ABOVE = 15
EMPTY_MSG = "Nenhum candidato MIN caiu dentro de internal_dark"


class StageError(RuntimeError):
    pass


def one(directory: Path, pattern: str) -> Path:
    hits = sorted(Path(directory).glob(pattern))
    if len(hits) != 1:
        raise RuntimeError(f"expected one {pattern} in {directory}, found {len(hits)}")
    return hits[0]


class Runner:
    def __init__(self, out: Path):
        self.out = out
        self.log = out / "commands.log"
        self.steps: list[dict] = []

    def __call__(self, name: str, cmd: list) -> int:
        cmd = [str(c) for c in cmd]
        t0 = time.time()
        rc, peak = run_tracked(cmd, self.log)
        self.steps.append({"step": name, "cmd0": cmd[2] if len(cmd) > 2 else cmd[0], "seconds": round(time.time() - t0, 1),
                           "rc": rc, "peak_memory_mb": round(peak, 1)})
        print(f"[v3] {self.out.name} {name}: rc={rc} {time.time() - t0:.0f}s peak {peak:.0f} MB", flush=True)
        return rc

    def script(self, name: str, stage: str, args: list, ok_empty: bool = False) -> str:
        rc = self(name, [PY, "-I", SCRIPT[stage], *args])
        if rc == 0:
            return "ok"
        tail = self.log.read_text(encoding="utf-8", errors="replace")[-4000:]
        if ok_empty and EMPTY_MSG in tail:
            return "empty"
        raise StageError(f"{name} failed rc={rc}; see {self.log}")


def run_crop(crop_id: str, gray: Path, out: Path, preset: str, e13_params: Path | None = None) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    run = Runner(out)
    rec: dict = {"crop_id": crop_id, "preset": preset, "gray": str(gray), "git": git_state(),
                 "timestamp": dt.datetime.now().astimezone().isoformat(timespec="seconds")}
    # 0 scale
    specimen = None
    if preset == "v3":
        sdir = out / "00_scale"
        sdir.mkdir(exist_ok=True)
        rc = run("scale", [PY, "-I", "-m", "batos_v3.scale", gray, "--out-json", sdir / "scale.json",
                           "--specimen-tif", sdir / "specimen.tif"])
        if rc:
            raise StageError("scale failed")
        sc = json.loads((sdir / "scale.json").read_text(encoding="utf-8"))
        r_scale, radius = float(sc["r_granulo"]), int(sc["openrec_radius"])
        s = r_scale / constants.R_REF
        flags = constants.stage_flags(s)
        specimen = tifffile.imread(sdir / "specimen.tif").astype(bool)
        rec |= {"r_granulo": r_scale, "r0": sc["r0"], "specimen_fraction": sc["specimen_fraction"], "s": s,
                "size_filter": constants.size_filter(s)}
    elif preset == "fixflat":
        params = json.loads(Path(e13_params).read_text(encoding="utf-8"))
        r_scale, radius, s = float(params["estimated_grain_radius"]), int(params["openrec_radius"]), None
        flags = {k: [str(x) for x in v] for k, v in BASE_FLAGS.items()}
        flags["flatten"] = flags["flatten"] + ["--variant", "A"]
        rec |= {"r_est_e13": r_scale, "e13_params": str(e13_params), "size_filter": 500}
    else:
        raise ValueError(preset)
    rec |= {"openrec_radius": radius, "flags": flags}
    # 1 openrec
    impl = "openrec_lowmem" if radius > LOWMEM_ABOVE else "openrec"
    prep = out / "01_openrec"
    run.script("openrec", impl, [gray, "--out-dir", prep, "--prefix", crop_id, "--radius", radius])
    flat = one(prep, "*_flatten_openrec_r*_8u.tif")
    rec["openrec_implementation"] = impl
    # 2 flatten
    pipe = out / "pipe"
    stem = flat.stem
    if run("flatten", [PY, "-I", "-m", "batos_v2.flatten_fixed", flat, "--out-dir", pipe / "01_preprocess",
                       "--prefix", stem, *flags["flatten"]]):
        raise StageError("flatten failed")
    selector = one(pipe / "01_preprocess", f"{stem}_flatten_lambda*_8u.tif")
    rec["bg_ref"] = json.loads((pipe / "01_preprocess" / "flatten_fixed_info.json").read_text(encoding="utf-8"))["bg_ref"]
    # 3 legacy ToS seed
    run.script("tos_seed", "tos_seed", [selector, "--out-dir", pipe / "02_tos_seed", *flags["tos_seed"]])
    raw_seed = one(pipe / "02_tos_seed", "*_tosconop_autodual_seed_raw.tif")
    # 4 background
    bg_flags = list(flags["background"])
    if preset == "v3":
        sel = tifffile.imread(selector)
        thr = otsu_u8(sel[specimen])
        bg_flags += ["--dark-threshold", str(thr)]
        rec["stage03_dark_threshold"] = thr
    run.script("background", "background", [flat, "--flattened", selector, "--seed", raw_seed,
                                            "--out-dir", pipe / "03_background", "--prefix", stem, *bg_flags])
    internal = one(pipe / "03_background", "*_internal_dark_candidates.tif")
    # 5 legacy markers
    run.script("markers_legacy", "markers_legacy", [flat, "--flattened", selector, "--internal-dark", internal,
                                                    "--out-dir", pipe / "04_markers", "--prefix", stem,
                                                    *flags["markers_legacy"]])
    baseline = one(pipe / "04_markers", "*_batos_marker_components_raw.tif")
    # 6 true dual, 7 spatialize
    e11 = out / "e11"
    run.script("true_dual", "true_dual", [selector, "--out-dir", e11 / "selector", *flags["true_dual"]])
    cand = one(e11 / "selector", "*_tosconop_autodual_candidates.csv")
    run.script("spatialize", "spatialize", ["--image", selector, "--candidate-csv", cand, "--internal-source", internal,
                                            "--out-dir", e11 / "spatial"])
    diag = e11 / "spatial" / "candidate_diagnostics_no_gt.csv"
    # 8 ToS-MIN markers
    status = run.script("markers_tos", "markers_tos", ["--diagnostic", diag, "--internal-source", internal,
                                                       *flags["markers_tos"], "--out-dir", e11 / "markers_tos_min"],
                        ok_empty=(preset == "v3"))
    if status == "empty":
        rec |= {"empty_result": f"e11 markers_tos: {EMPTY_MSG} (no MIN candidate inside internal_dark)",
                "steps": run.steps, "seconds_total": round(time.time() - t_start, 1)}
        rec["peak_memory_mb"] = max(st["peak_memory_mb"] for st in run.steps)
        (out / "v3_run.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
        (out / "EMPTY").write_text(rec["empty_result"], encoding="utf-8")
        return rec
    tos_dir = e11 / "markers_tos_min"
    tos_marker = tos_dir / "markers_uint16.tif"
    # 9 EDT-Otsu control
    run.script("edt_otsu", "edt_otsu", ["--image", flat, *flags["edt_otsu"], "--match-count-from", tos_marker,
                                        "--out-dir", e11 / "markers_edt_otsu"])
    # 10 watershed, ToS route
    run.script("ws_tos", "ws_tos", ["--image", flat, "--mask-marker", tos_marker, "--mask-marker", baseline,
                                    *flags["ws_tos"], "--marker", f"tos_e11_min={tos_marker}",
                                    "--marker", f"baseline_existing={baseline}", "--out-dir", e11 / "watershed"])
    # 11 E23c route, 12 watershed selective
    e23 = out / "e23c"
    run.script("route", "route", ["--candidate-diagnostics", diag, "--tos-marker-csv", tos_dir / "markers.csv",
                                  "--tos-marker-tiff", tos_marker,
                                  "--edt-marker-tiff", e11 / "markers_edt_otsu" / "markers_uint16.tif",
                                  "--crop-id", crop_id, "--out-dir", e23 / "route", *flags.get("route", [])])
    run.script("ws_sel", "ws_sel", ["--image", flat, "--mask-marker", tos_marker,
                                    "--mask-marker", e11 / "watershed" / "normalized_markers" / "baseline_existing.tif",
                                    *flags["ws_sel"],
                                    "--marker", f"batos_selective={e23 / 'route' / 'batos_selective_markers_uint16.tif'}",
                                    "--out-dir", e23 / "watershed"])
    gate = json.loads((e23 / "route" / "route_manifest.json").read_text(encoding="utf-8"))
    gate = gate.get("gate", gate)
    rec |= {"gate": {k: gate.get(k) for k in ("gate_triggered", "region_voxels_median", "region_voxels_q75",
                                               "median_limit", "q75_limit")},
            "tos_markers": int(sum(1 for _ in open(tos_dir / "markers.csv", encoding="utf-8")) - 1),
            "steps": run.steps, "seconds_total": round(time.time() - t_start, 1),
            "peak_memory_mb": max(st["peak_memory_mb"] for st in run.steps)}
    (out / "v3_run.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    (out / "DONE").write_text("ok", encoding="utf-8")
    return rec


SEGMENTATIONS = {"BA-TOS_ToS": ("e11/watershed", "tos_e11_min"), "BA-TOS_sel": ("e23c/watershed", "batos_selective")}


def outputs(out: Path, shape) -> dict:
    """Per route (segmentation, normalized markers); empty arrays for an EMPTY result."""
    res = {}
    for route, (ws, label) in SEGMENTATIONS.items():
        if (out / "EMPTY").exists():
            z = np.zeros(shape, np.uint32)
            res[route] = (z, z)
        else:
            res[route] = (tifffile.imread(out / ws / "segmentations" / f"{label}_distance_c1.tif"),
                          tifffile.imread(out / ws / "normalized_markers" / f"{label}.tif"))
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--crop-id", required=True)
    ap.add_argument("--gray", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--preset", choices=["v3", "fixflat"], required=True)
    ap.add_argument("--e13-params", type=Path)
    a = ap.parse_args(argv)
    rec = run_crop(a.crop_id, a.gray, a.out, a.preset, a.e13_params)
    print(json.dumps({k: rec.get(k) for k in ("crop_id", "preset", "openrec_radius", "s", "stage03_dark_threshold",
                                              "gate", "tos_markers", "empty_result", "seconds_total")}, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
