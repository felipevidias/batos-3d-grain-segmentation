"""Phase 1c round 4: run the FROZEN E-series code verbatim on one crop.

Every stage is an unmodified recovered script called as a subprocess (python -I), with the flags
and stage order of e13_stage2_worker_20260725_165544.sh -> run_e13_e11_all_crops.sh (chat_recovered)
and of the E13 preprocessing (e13_reconstruct_official_pipeline.py / run_batos_openrec_auto_minparam.py).
The SHA-256 of every script is recorded in frozen_run.json. Nothing is re-implemented here; this
module only builds command lines and checks that expected outputs exist (fail loudly).

Stage order (one crop):
  0 scale      e13_estimate_scale.py (12-row manifest, GT-free)    -> R_hat, openrec_radius, min_size
  1 openrec    batos_preprocess_opening_reconstruction_v2.py --radius r   (zsh_recovered)
  2 pipeline   run_batos_pipeline.sh steps 1-4 with lambda 1 (python steps called individually,
               flags identical to the recovered script): batos_preprocess.py, batos_tos_seed.py,
               batos_background.py, batos_markers.py. Steps 5-6 (legacy watershed) are not consumed by E11.
  3 selector   batos_tos_seed_true_dual.py (flags of run_e13_e11_all_crops.sh)
  4 spatial    e10_spatialize_candidates_no_gt.py
  5 markers    e11_build_markers_ablation.py (min, dual), e11_build_control_markers.py (otsu_edt)
  6 watershed  e11_run_fixed_mask_watershed.py (mask markers: ToS MIN + baseline)
  7 evaluate   e09b_ignore_aware_audit.py
  8 route      e23c_build_selective_marker_route.py -> watershed + e09b for batos_selective (E23c runner)
Not run: audit_true_batos_e07.py (GT capacity audit; not recovered), internal_medial / mintree_e14 controls.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from .common import REPO_ROOT
from .run_crop import git_state

LI = Path(r"C:\BA_TOS_DATA\linux_import")
CHAT = LI / "chat_recovered" / "e13_crop05_metadata_full_20260725_182408" / "audit"
SCRIPTS = {
    "e13_estimate_scale": LI / "linux_patches/BA_TOS_E13_E16_conditional_acceptance/audit/e13_estimate_scale.py",
    "openrec_v2": LI / "zsh_recovered/code__segmentation__batos__batos_preprocess_opening_reconstruction_v2.py",
    "batos_preprocess": REPO_ROOT / "src/batos/batos_preprocess.py",
    "batos_tos_seed": Path(r"C:\BA_TOS_DATA\frozen_code\pilot_branch_413920f\batos_tos_seed.py"),
    "batos_background": REPO_ROOT / "src/batos/batos_background.py",
    "batos_markers": REPO_ROOT / "src/batos/batos_markers.py",
    "true_dual": CHAT / "batos_tos_seed_true_dual.py",
    "spatialize": CHAT / "e10_spatialize_candidates_no_gt.py",
    "e11_markers": CHAT / "e11_build_markers_ablation.py",
    "e11_controls": CHAT / "e11_build_control_markers.py",
    "e11_watershed": CHAT / "e11_run_fixed_mask_watershed.py",
    "e09b": CHAT / "e09b_ignore_aware_audit.py",
    "e23c_route": LI / "linux_patches/e23c_selective_router_bundle/e23c_build_selective_marker_route.py",
}

# run_batos_pipeline.sh (recovered) step flags
TOS_SEED_LEGACY = ["--height-factor", "0.012", "--area-factor", "22", "--area-min", "2", "--area-max", "2200",
                   "--delta-parent-min", "0", "--dark-percentile-max", "78", "--min-score", "0.35",
                   "--max-selected-nodes", "900", "--support-area-min", "20", "--support-area-max", "45000",
                   "--support-max-steps", "16", "--component-min", "2", "--component-max", "1500"]
BACKGROUND = ["--dark-percentile", "40", "--otsu-offset", "-10", "--dark-threshold-max", "135",
              "--dark-connectivity", "1", "--seed-area-min", "1", "--seed-area-max", "3000",
              "--max-external-overlap", "0.10", "--min-internal-overlap", "0.00", "--internal-search-radius", "6",
              "--center-dark-max", "210", "--min-ring-contrast", "-5", "--max-ring-external-coverage", "0.70",
              "--nms-min-dist", "5"]
MARKERS_LEGACY = ["--connectivity", "1", "--area-min", "3", "--area-max", "900", "--max-dim", "30",
                  "--max-elongation", "8", "--min-fill", "0.005", "--center-dark-max", "205",
                  "--min-ring-contrast", "-2", "--nms-min-dist", "7", "--max-markers", "420"]
# run_e13_e11_all_crops.sh selector flags
TRUE_DUAL = ["--height-factor", "0.012", "--area-factor", "22.0", "--area-min", "2", "--area-max", "2200",
             "--delta-parent-min", "0", "--dark-percentile-max", "78.0", "--support-area-min", "20",
             "--support-area-max", "45000", "--support-max-steps", "16", "--support-area-ratio-preferred", "18.0",
             "--support-area-ratio-width", "80.0", "--support-contrast-norm", "80.0", "--delta-norm", "20.0",
             "--area-preferred-position", "0.42", "--area-preferred-width", "0.55", "--w-height", "0.22",
             "--w-delta", "0.12", "--w-area", "0.18", "--w-dark", "0.18", "--w-support", "0.3",
             "--min-score", "0.35", "--cleanup-open-radius", "0", "--cleanup-close-radius", "0",
             "--component-min", "2", "--component-max", "1500", "--max-selected-nodes", "3600"]
MARKER_FLAGS = ["--nms-radius", "6", "--edt-priority-weight", "0.20", "--min-priority-bonus", "0.05"]


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def one(directory: Path, pattern: str) -> Path:
    hits = sorted(Path(directory).glob(pattern))
    if len(hits) != 1:
        raise RuntimeError(f"expected one {pattern} in {directory}, found {len(hits)}: {hits}")
    return hits[0]


class Runner:
    def __init__(self, log_path: Path):
        self.log_path = log_path
        self.steps: list[dict] = []

    def __call__(self, name: str, script_key: str, args: list[str], cwd: Path | None = None) -> None:
        script = SCRIPTS[script_key]
        if not script.is_file():
            raise FileNotFoundError(script)
        cmd = [sys.executable, "-I", str(script), *[str(a) for a in args]]
        t0 = time.time()
        with open(self.log_path, "a", encoding="utf-8") as log:
            log.write(f"\n$ {subprocess.list2cmdline(cmd)}\n")
            log.flush()
            import psutil
            proc = psutil.Popen(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT)
            peak = 0
            while proc.poll() is None:  # peak working set of the stage subprocess (Windows: peak_wset)
                try:
                    info = proc.memory_info()
                    peak = max(peak, getattr(info, "peak_wset", 0) or info.rss)
                except psutil.Error:
                    pass
                time.sleep(0.2)
            returncode = proc.returncode
        self.steps.append({"step": name, "script": str(script), "sha256": sha256(script),
                           "command": subprocess.list2cmdline(cmd), "seconds": round(time.time() - t0, 1),
                           "peak_memory_mb": round(peak / 2 ** 20, 1), "returncode": returncode})
        print(f"[frozen] {name}: rc={returncode} ({time.time() - t0:.0f}s, peak {peak / 2 ** 20:.0f} MB)", flush=True)
        if returncode != 0:
            raise RuntimeError(f"{name} failed; see {self.log_path}")


def scale_params(crops_meta: dict, out_dir: Path, run: Runner) -> Path:
    """e13_estimate_scale.py requires exactly 12 manifest rows (scan, crop_id, gray)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = out_dir / "scale_manifest.csv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scan", "crop_id", "gray"])
        writer.writeheader()
        for crop in crops_meta["crops"]:
            writer.writerow({"scan": crop["scan"], "crop_id": crop["crop_id"],
                             "gray": crop["gray_crop_path"]})
    params = out_dir / "auto_scale_params"
    if not (params / "e13_auto_scale_params_per_crop.csv").exists():
        run("scale", "e13_estimate_scale", ["--manifest", manifest, "--out-dir", params])
    return params


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Run the frozen E-series code verbatim on one crop")
    parser.add_argument("--crop-id", required=True)
    parser.add_argument("--runs-root", type=Path, default=Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen"))
    parser.add_argument("--skip-route", action="store_true")
    args = parser.parse_args(argv)

    meta = json.loads((REPO_ROOT / "data/metadata/crops_sibgrapi.json").read_text(encoding="utf-8"))
    crop = next(c for c in meta["crops"] if c["crop_id"] == args.crop_id)
    root = args.runs_root
    out = root / args.crop_id
    if (out / "frozen_run.json").exists():
        raise FileExistsError(f"{out} already completed; use a new --runs-root")
    out.mkdir(parents=True, exist_ok=True)
    run = Runner(out / "frozen_commands.log")
    t_start = time.time()
    cid = args.crop_id
    gray, gt = Path(crop["gray_crop_path"]), Path(crop["gt_path"])

    # 0 scale (shared across crops)
    params_dir = scale_params(meta, root / "_scale", run)
    params = json.loads((params_dir / f"{cid}_auto_params.json").read_text(encoding="utf-8"))
    radius, min_size = int(params["openrec_radius"]), int(params["min_predicted_size"])

    # 1 openrec v2 (E13: --prefix crop_id)
    prep = out / "01_preprocess_openrec_auto_minparam"
    run("openrec_v2", "openrec_v2", [gray, "--out-dir", prep, "--prefix", cid, "--radius", radius])
    flat_image = one(prep, "*_flatten_openrec_r*_8u.tif")

    # 2 run_batos_pipeline.sh <flat> <out> 1 <min_size>: steps 1-4
    pipe = out / "batos_openrec_auto_minparam"
    stem = flat_image.stem
    run("pipeline_1_flatten_lambda1", "batos_preprocess",
        [flat_image, "--out-dir", pipe / "01_preprocess", "--prefix", stem, "--lambda-xy", 1])
    selector = one(pipe / "01_preprocess", f"{stem}_flatten_lambda*_8u.tif")
    run("pipeline_2_tos_seed_k900", "batos_tos_seed", [selector, "--out-dir", pipe / "02_tos_seed", *TOS_SEED_LEGACY])
    raw_seed = one(pipe / "02_tos_seed", "*_tosconop_autodual_seed_raw.tif")
    run("pipeline_3_background", "batos_background",
        [flat_image, "--flattened", selector, "--seed", raw_seed, "--out-dir", pipe / "03_background",
         "--prefix", stem, *BACKGROUND])
    internal = one(pipe / "03_background", "*_internal_dark_candidates.tif")
    run("pipeline_4_markers", "batos_markers",
        [flat_image, "--flattened", selector, "--internal-dark", internal, "--out-dir", pipe / "04_markers",
         "--prefix", stem, *MARKERS_LEGACY])
    baseline = one(pipe / "04_markers", "*_batos_marker_components_raw.tif")

    # 3-7 E11 (run_e13_e11_all_crops.sh order)
    e11 = out / "e11"
    run("selector_true_dual_k3600", "true_dual", [selector, "--out-dir", e11 / "selector_k3600", *TRUE_DUAL])
    cand_csv = one(e11 / "selector_k3600", "*_tosconop_autodual_candidates.csv")
    run("spatialize", "spatialize", ["--image", selector, "--candidate-csv", cand_csv,
                                     "--internal-source", internal, "--out-dir", e11 / "spatial_no_gt"])
    diag = e11 / "spatial_no_gt" / "candidate_diagnostics_no_gt.csv"
    run("markers_tos_min", "e11_markers", ["--diagnostic", diag, "--internal-source", internal,
                                           "--polarity", "min", *MARKER_FLAGS, "--out-dir", e11 / "markers_tos_min"])
    tos_marker = e11 / "markers_tos_min" / "markers_uint16.tif"
    run("markers_tos_dual", "e11_markers", ["--diagnostic", diag, "--internal-source", internal,
                                            "--polarity", "dual", *MARKER_FLAGS, "--out-dir", e11 / "markers_tos_dual"])
    run("markers_edt_otsu", "e11_controls", ["--image", flat_image, "--mode", "otsu_edt", "--nms-radius", "6",
                                             "--match-count-from", tos_marker, "--out-dir", e11 / "markers_edt_otsu"])
    ws_args = ["--image", flat_image, "--mask-marker", tos_marker, "--mask-marker", baseline,
               "--outside-policy", "drop", "--marker", f"tos_e11_min={tos_marker}",
               "--marker", f"baseline_existing={baseline}",
               "--marker", f"edt_otsu={e11 / 'markers_edt_otsu' / 'markers_uint16.tif'}",
               "--marker", f"tos_e11_dual={e11 / 'markers_tos_dual' / 'markers_uint16.tif'}",
               "--out-dir", e11 / "watershed"]
    run("watershed_fixed_mask", "e11_watershed", ws_args)
    ws = pd.read_csv(e11 / "watershed" / "watershed_results.csv")
    eval_args = ["--e09-dir", e11 / "watershed", "--gt", gt, "--tos-label", "tos_e11_min",
                 "--out-dir", e11 / "evaluation"]
    for _, r in ws.iterrows():
        eval_args += ["--marker", f"{r['marker_source']}={r['normalized_marker']}"]
    run("evaluate_e09b", "e09b", eval_args)

    # 8 E23c route
    route = {}
    if not args.skip_route:
        e23 = out / "e23c"
        run("route_e23c", "e23c_route", ["--candidate-diagnostics", diag,
                                         "--tos-marker-csv", e11 / "markers_tos_min" / "markers.csv",
                                         "--tos-marker-tiff", tos_marker,
                                         "--edt-marker-tiff", e11 / "markers_edt_otsu" / "markers_uint16.tif",
                                         "--crop-id", cid, "--out-dir", e23 / "route"])
        baseline_norm = e11 / "watershed" / "normalized_markers" / "baseline_existing.tif"
        run("watershed_selective", "e11_watershed",
            ["--image", flat_image, "--mask-marker", tos_marker, "--mask-marker", baseline_norm,
             "--outside-policy", "drop",
             "--marker", f"batos_selective={e23 / 'route' / 'batos_selective_markers_uint16.tif'}",
             "--out-dir", e23 / "watershed"])
        run("evaluate_selective_e09b", "e09b",
            ["--e09-dir", e23 / "watershed", "--gt", gt,
             "--marker", f"batos_selective={e23 / 'watershed' / 'normalized_markers' / 'batos_selective.tif'}",
             "--tos-label", "batos_selective", "--out-dir", e23 / "evaluation"])
        route = json.loads((e23 / "route" / "route_manifest.json").read_text(encoding="utf-8"))

    results = pd.read_csv(e11 / "evaluation" / "ignore_aware_results.csv")
    record = {
        "experiment": "batos-paper-v2-frozen", "crop_id": cid,
        "timestamp": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "command": subprocess.list2cmdline([sys.executable, "-I", "-m", "batos_v2.frozen", *(argv or sys.argv[1:])]),
        "git": git_state(), "auto_params": params, "openrec_radius": radius, "min_predicted_size": min_size,
        "inputs": {"gray": str(gray), "gt": str(gt)},
        "outputs": {"flat_image": str(flat_image), "selector": str(selector), "internal": str(internal),
                    "baseline_markers": str(baseline), "candidates": str(cand_csv), "diagnostics": str(diag)},
        "steps": run.steps, "route": route,
        "e09b_results": results.to_dict(orient="records"),
        "seconds_total": round(time.time() - t_start, 1),
        "peak_memory_mb": max((s.get("peak_memory_mb", 0) for s in run.steps), default=None),
    }
    (out / "frozen_run.json").write_text(json.dumps(record, indent=2, default=str), encoding="utf-8")
    cols = [c for c in ("evaluation_protocol", "marker_source", "f1_iou025", "f1_iou050", "f1_iou075",
                        "precision_iou050", "recall_iou050") if c in results.columns]
    print(results[cols].to_string(index=False))
    if route:
        g = route.get("gate", route)
        print(f"[frozen] gate: {json.dumps(g)[:300]}")
    print(f"[frozen] done in {record['seconds_total']:.0f}s -> {out / 'frozen_run.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
