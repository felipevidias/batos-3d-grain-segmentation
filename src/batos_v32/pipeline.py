"""BA-TOS v3.2 single-crop driver (Phase 5h): v3.1 stages up to the e11 ToS-MIN markers (prefix, run once per crop),
then a marker-uniqueness rule (batos_v32.uniqueness) and the unchanged v3 tail (EDT-Otsu control matched to the new
count, fixed-mask watershed, E23c route, selective watershed), once per hypothesis.

The stage calls are those of batos_v3.pipeline.run_crop (preset v3) with the v3.1 count rule (batos_v31.constants),
split in two so that the expensive prefix is shared by all hypotheses. With hypothesis `none` the tail reproduces v3.1
(checked bit-for-bit). No v3 / v3.1 file is modified. GT is never read.
"""

from __future__ import annotations

import datetime as dt
import json
import sys
import time
from pathlib import Path

import tifffile

from batos_v2.run_crop import git_state
from batos_v3.pipeline import EMPTY_MSG, PY, Runner, StageError, one
from batos_v3.scale import otsu_u8
from batos_v31 import constants as C31

HERE = Path(__file__).resolve().parent
PREFIX_DONE = "PREFIX_DONE"


def run_prefix(crop_id: str, gray: Path, out: Path) -> dict:
    """Stages 0-8 of batos_v3.pipeline.run_crop (preset v3) with v3.1 flags. Writes prefix.json."""
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    run = Runner(out)
    rec: dict = {"crop_id": crop_id, "preset": "v3.1-prefix", "gray": str(gray), "git": git_state(),
                 "timestamp": dt.datetime.now().astimezone().isoformat(timespec="seconds")}
    sdir = out / "00_scale"
    sdir.mkdir(exist_ok=True)
    if run("scale", [PY, "-I", "-m", "batos_v3.scale", gray, "--out-json", sdir / "scale.json",
                     "--specimen-tif", sdir / "specimen.tif"]):
        raise StageError("scale failed")
    sc = json.loads((sdir / "scale.json").read_text(encoding="utf-8"))
    r_scale, radius = float(sc["r_granulo"]), int(sc["openrec_radius"])
    s = r_scale / C31.R_REF
    flags = C31.stage_flags(s)
    specimen = tifffile.imread(sdir / "specimen.tif").astype(bool)
    rec |= {"r_granulo": r_scale, "r0": sc["r0"], "specimen_fraction": sc["specimen_fraction"], "s": s,
            "size_filter": C31.size_filter(s), "openrec_radius": radius, "flags": flags,
            "n_rule": "v2 * max(1, s^-3)", "n_values": {f"{c.stage} {c.flag}": C31.scaled(c, s) for c in C31.N_CONSTANTS}}
    impl = "openrec_lowmem" if radius > 15 else "openrec"
    prep = out / "01_openrec"
    run.script("openrec", impl, [gray, "--out-dir", prep, "--prefix", crop_id, "--radius", radius])
    flat = one(prep, "*_flatten_openrec_r*_8u.tif")
    rec["openrec_implementation"] = impl
    pipe = out / "pipe"
    stem = flat.stem
    if run("flatten", [PY, "-I", "-m", "batos_v2.flatten_fixed", flat, "--out-dir", pipe / "01_preprocess",
                       "--prefix", stem, *flags["flatten"]]):
        raise StageError("flatten failed")
    selector = one(pipe / "01_preprocess", f"{stem}_flatten_lambda*_8u.tif")
    rec["bg_ref"] = json.loads((pipe / "01_preprocess" / "flatten_fixed_info.json").read_text(encoding="utf-8"))["bg_ref"]
    run.script("tos_seed", "tos_seed", [selector, "--out-dir", pipe / "02_tos_seed", *flags["tos_seed"]])
    raw_seed = one(pipe / "02_tos_seed", "*_tosconop_autodual_seed_raw.tif")
    thr = otsu_u8(tifffile.imread(selector)[specimen])
    rec["stage03_dark_threshold"] = thr
    run.script("background", "background", [flat, "--flattened", selector, "--seed", raw_seed,
                                            "--out-dir", pipe / "03_background", "--prefix", stem,
                                            *flags["background"], "--dark-threshold", str(thr)])
    internal = one(pipe / "03_background", "*_internal_dark_candidates.tif")
    run.script("markers_legacy", "markers_legacy", [flat, "--flattened", selector, "--internal-dark", internal,
                                                    "--out-dir", pipe / "04_markers", "--prefix", stem,
                                                    *flags["markers_legacy"]])
    baseline = one(pipe / "04_markers", "*_batos_marker_components_raw.tif")
    e11 = out / "e11"
    run.script("true_dual", "true_dual", [selector, "--out-dir", e11 / "selector", *flags["true_dual"]])
    cand = one(e11 / "selector", "*_tosconop_autodual_candidates.csv")
    run.script("spatialize", "spatialize", ["--image", selector, "--candidate-csv", cand, "--internal-source", internal,
                                            "--out-dir", e11 / "spatial"])
    status = run.script("markers_tos", "markers_tos", ["--diagnostic", e11 / "spatial" / "candidate_diagnostics_no_gt.csv",
                                                       "--internal-source", internal, *flags["markers_tos"],
                                                       "--out-dir", e11 / "markers_tos_min"], ok_empty=True)
    rec |= {"paths": {"flat": flat, "selector": selector, "internal": internal, "baseline": baseline, "candidates": cand,
                      "diag": e11 / "spatial" / "candidate_diagnostics_no_gt.csv", "tos_dir": e11 / "markers_tos_min"},
            "empty_result": f"e11 markers_tos: {EMPTY_MSG}" if status == "empty" else None,
            "steps": run.steps, "seconds_prefix": round(time.time() - t0, 1),
            "peak_memory_mb": max(st["peak_memory_mb"] for st in run.steps)}
    (out / "prefix.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    (out / PREFIX_DONE).write_text("ok", encoding="utf-8")
    return rec


def ensure_supports(prefix: Path) -> Path:
    """Exclusive ToS supports (H3 input), computed once per crop in a subprocess (peak ~3 GB)."""
    p = prefix / "e11" / "supports" / "leaf_labels.tif"
    if p.exists():
        return p
    rec = json.loads((prefix / "prefix.json").read_text(encoding="utf-8"))
    run = Runner(prefix)
    if run("supports", [PY, "-I", HERE / "supports.py", "--image", rec["paths"]["selector"],
                        "--candidate-csv", rec["paths"]["candidates"], "--out", p]):
        raise StageError("supports failed")
    return p


def run_tail(prefix: Path, hyp: str, out: Path) -> dict:
    """Uniqueness rule + stages 9-12 of batos_v3.pipeline.run_crop. Writes v3_run.json (+ DONE / EMPTY)."""
    from batos_v32 import uniqueness
    pre = json.loads((prefix / "prefix.json").read_text(encoding="utf-8"))
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    run = Runner(out)
    flags, s, r = pre["flags"], pre["s"], pre["r_granulo"]
    P = {k: Path(v) for k, v in pre["paths"].items()}
    rec = {k: pre[k] for k in ("crop_id", "gray", "r_granulo", "r0", "specimen_fraction", "s", "size_filter",
                               "openrec_radius", "flags", "n_rule", "n_values", "stage03_dark_threshold", "bg_ref",
                               "openrec_implementation")}
    rec |= {"preset": f"v3.2-{hyp}", "hypothesis": hyp, "prefix": str(prefix), "git": git_state(),
            "timestamp": dt.datetime.now().astimezone().isoformat(timespec="seconds")}
    if pre.get("empty_result"):
        rec |= {"empty_result": pre["empty_result"], "steps": pre["steps"], "seconds_total": pre["seconds_prefix"],
                "peak_memory_mb": pre["peak_memory_mb"]}
        (out / "v3_run.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
        (out / "EMPTY").write_text(rec["empty_result"], encoding="utf-8")
        return rec
    e11 = out / "e11"
    tos_dir = e11 / "markers_tos_unique"
    fl = dict(zip(flags["edt_otsu"][::2], flags["edt_otsu"][1::2]))
    kw = {"flat": P["flat"], "diag": P["diag"], "min_component": int(fl["--v3-min-component"]),
          "closing": int(fl["--v3-closing-iterations"]), "edt_nms": float(fl["--nms-radius"])}
    if hyp == "H3":
        kw["leaf"] = ensure_supports(prefix)
    tu = time.time()
    rec["uniqueness"] = uniqueness.apply(hyp, P["tos_dir"], tos_dir, r_granulo=r, s=s, **kw)
    run.steps.append({"step": "uniqueness", "seconds": round(time.time() - tu, 1), "rc": 0})
    tos_marker = tos_dir / "markers_uint16.tif"
    if rec["uniqueness"]["markers_out"] == 0:
        raise StageError("uniqueness removed every marker")
    run.script("edt_otsu", "edt_otsu", ["--image", P["flat"], *flags["edt_otsu"], "--match-count-from", tos_marker,
                                        "--out-dir", e11 / "markers_edt_otsu"])
    run.script("ws_tos", "ws_tos", ["--image", P["flat"], "--mask-marker", tos_marker, "--mask-marker", P["baseline"],
                                    *flags["ws_tos"], "--marker", f"tos_e11_min={tos_marker}",
                                    "--marker", f"baseline_existing={P['baseline']}", "--out-dir", e11 / "watershed"])
    e23 = out / "e23c"
    run.script("route", "route", ["--candidate-diagnostics", P["diag"], "--tos-marker-csv", tos_dir / "markers.csv",
                                  "--tos-marker-tiff", tos_marker,
                                  "--edt-marker-tiff", e11 / "markers_edt_otsu" / "markers_uint16.tif",
                                  "--crop-id", pre["crop_id"], "--out-dir", e23 / "route", *flags.get("route", [])])
    run.script("ws_sel", "ws_sel", ["--image", P["flat"], "--mask-marker", tos_marker,
                                    "--mask-marker", e11 / "watershed" / "normalized_markers" / "baseline_existing.tif",
                                    *flags["ws_sel"],
                                    "--marker", f"batos_selective={e23 / 'route' / 'batos_selective_markers_uint16.tif'}",
                                    "--out-dir", e23 / "watershed"])
    gate = json.loads((e23 / "route" / "route_manifest.json").read_text(encoding="utf-8"))
    gate = gate.get("gate", gate)
    rec |= {"gate": {k: gate.get(k) for k in ("gate_triggered", "region_voxels_median", "region_voxels_q75",
                                               "median_limit", "q75_limit")},
            "tos_markers_prefix": rec["uniqueness"]["markers_in"], "tos_markers": rec["uniqueness"]["markers_out"],
            "steps": pre["steps"] + run.steps, "seconds_tail": round(time.time() - t0, 1),
            "seconds_total": round(pre["seconds_prefix"] + time.time() - t0, 1),
            "peak_memory_mb": max(st.get("peak_memory_mb", 0) for st in pre["steps"] + run.steps)}
    (out / "v3_run.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    (out / "DONE").write_text("ok", encoding="utf-8")
    return rec


if __name__ == "__main__":
    sys.exit("use scripts/run_v32_dev.py")
