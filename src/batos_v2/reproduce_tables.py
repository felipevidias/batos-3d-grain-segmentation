"""Phase 1e: reproduce paper Tables II (Min-tree baseline), IV (sensitivity) and V (crop extent) with the frozen code.

All stages are unmodified recovered scripts, run through batos_v2.frozen.Runner (subprocess, python -I,
SHA-256 + runtime + process-tree peak memory per step). The frozen 12-crop runs in
C:\\BA_TOS_DATA\\runs\\batos_paper_v2_frozen are read-only inputs; every output goes to --root
(default C:\\BA_TOS_DATA\\runs\\reproduction_tables) and a subcommand refuses to overwrite a finished run.

  mintree --crop-id ID   E14 Min-tree source control, flags of run_e13_e11_all_crops.sh:
                         e14_build_mintree_candidates.py (--max-candidates 3600 --area-min 20 --area-max 45000)
                         -> e11 markers (min, NMS 6) -> fixed-mask watershed (mask: ToS MIN + baseline, as in the
                         original single call) with markers tos_e11_min (consistency check) + mintree_e14 -> e09b.
  sensitivity            run_e13_development_sensitivity.sh on EFRGP01_00_crop04: tau_s 0.25..0.45 (K3600, NMS6),
                         K 900/1800/7200 (tau 0.35, NMS6), NMS 4/8/10 (default spatial diagnostics); one fixed-mask
                         watershed (mask: default K3600/NMS6 markers + baseline) -> e09b.
  cropsize               run_e15_crop_size_sensitivity.sh: e15_make_crop_size_suite.py (centre crops 128/160/200 of
                         both crop04 parents, parent scale params frozen) -> batos_v2.frozen per sub-crop.
e09b is called through scripts/e09b_with_iou025.py (adds IoU 0.25; 0.50/0.75 identical to plain e09b).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from .common import REPO_ROOT
from .frozen import MARKER_FLAGS, SCRIPTS, TRUE_DUAL, Runner
from .run_crop import git_state

RUNS = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen")
DEFAULT_ROOT = Path(r"C:\BA_TOS_DATA\runs\reproduction_tables")
DEV_CROP = "EFRGP01_00_crop04"


def frozen_inputs(crop_id: str) -> dict:
    rec_path = RUNS / crop_id / "frozen_run.json"
    if not rec_path.is_file():
        raise FileNotFoundError(rec_path)
    rec = json.loads(rec_path.read_text(encoding="utf-8"))
    o = rec["outputs"]
    paths = {"selector": Path(o["selector"]), "internal": Path(o["internal"]), "flat": Path(o["flat_image"]),
             "baseline": Path(o["baseline_markers"]), "gt": Path(rec["inputs"]["gt"]),
             "tos_marker": RUNS / crop_id / "e11" / "markers_tos_min" / "markers_uint16.tif"}
    for key, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"{crop_id} {key}: {path}")
    return paths | {"frozen_run": rec_path}


def evaluate(run: Runner, name: str, ws_dir: Path, gt: Path, tos_label: str, out_dir: Path) -> pd.DataFrame:
    ws = pd.read_csv(ws_dir / "watershed_results.csv")
    args = [SCRIPTS["e09b"], "--e09-dir", ws_dir, "--gt", gt, "--tos-label", tos_label, "--out-dir", out_dir]
    for _, r in ws.iterrows():
        args += ["--marker", f"{r['marker_source']}={r['normalized_marker']}"]
    run(name, "e09b_iou025", args)
    return pd.read_csv(out_dir / "ignore_aware_results.csv")


def write_record(out: Path, experiment: str, run: Runner, t0: float, inputs: dict, results: pd.DataFrame,
                 extra: dict | None = None) -> None:
    record = {
        "experiment": experiment,
        "timestamp": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "command": subprocess.list2cmdline([sys.executable, "-I", "-m", "batos_v2.reproduce_tables", *sys.argv[1:]]),
        "git": git_state(), "inputs": {k: str(v) for k, v in inputs.items()},
        "steps": run.steps, "e09b_results": results.to_dict(orient="records"),
        "seconds_total": round(time.time() - t0, 1),
        "peak_memory_mb": max((s["peak_memory_mb"] for s in run.steps), default=None),
    } | (extra or {})
    (out / "run.json").write_text(json.dumps(record, indent=2, default=str), encoding="utf-8")


def start(out: Path) -> tuple[Runner, float]:
    if (out / "run.json").exists():
        raise FileExistsError(f"{out} already completed; use a new --root")
    out.mkdir(parents=True, exist_ok=True)
    return Runner(out / "commands.log"), time.time()


def mintree(crop_id: str, root: Path) -> None:
    src = frozen_inputs(crop_id)
    out = root / "mintree" / crop_id
    run, t0 = start(out)
    run("e14_mintree_candidates", "e14_mintree",
        ["--image", src["selector"], "--internal-source", src["internal"], "--max-candidates", 3600,
         "--area-min", 20, "--area-max", 45000, "--out-dir", out / "mintree_source_no_gt"])
    run("markers_mintree_e14", "e11_markers",
        ["--diagnostic", out / "mintree_source_no_gt" / "candidate_diagnostics_no_gt.csv",
         "--internal-source", src["internal"], "--polarity", "min", *MARKER_FLAGS,
         "--out-dir", out / "markers_mintree_e14"])
    run("watershed_fixed_mask", "e11_watershed",
        ["--image", src["flat"], "--mask-marker", src["tos_marker"], "--mask-marker", src["baseline"],
         "--outside-policy", "drop", "--marker", f"tos_e11_min={src['tos_marker']}",
         "--marker", f"mintree_e14={out / 'markers_mintree_e14' / 'markers_uint16.tif'}",
         "--out-dir", out / "watershed"])
    results = evaluate(run, "evaluate_e09b", out / "watershed", src["gt"], "tos_e11_min", out / "evaluation")
    write_record(out, "phase1e_table2_mintree_e14", run, t0, src, results, {"crop_id": crop_id})
    seed = results[results.evaluation_protocol == "seed_ignore"]
    print(seed[["marker_source", "f1_iou025", "f1_iou050", "f1_iou075"]].to_string(index=False), flush=True)


def sensitivity(root: Path) -> None:
    src = frozen_inputs(DEV_CROP)
    out = root / "sensitivity_crop04"
    run, t0 = start(out)
    variants = []

    def selector_variant(score: str, k: int, nms: int) -> None:
        label = f"score_{score.replace('.', 'p')}_k{k}_n{nms}"
        v = out / label
        flags = list(TRUE_DUAL)
        flags[flags.index("--min-score") + 1] = score
        flags[flags.index("--max-selected-nodes") + 1] = str(k)
        run(f"{label}:selector", "true_dual", [src["selector"], "--out-dir", v / "selector", *flags])
        hits = sorted((v / "selector").glob("*_tosconop_autodual_candidates.csv"))
        if len(hits) != 1:
            raise RuntimeError(f"{label}: expected one candidates CSV, found {hits}")
        run(f"{label}:spatialize", "spatialize", ["--image", src["selector"], "--candidate-csv", hits[0],
                                                  "--internal-source", src["internal"], "--out-dir", v / "spatial"])
        markers(label, v / "spatial" / "candidate_diagnostics_no_gt.csv", score, k, nms)

    def markers(label: str, diagnostic: Path, score: str, k: int, nms: int) -> None:
        v = out / label
        run(f"{label}:markers", "e11_markers",
            ["--diagnostic", diagnostic, "--internal-source", src["internal"], "--polarity", "min",
             "--nms-radius", nms, "--edt-priority-weight", "0.20", "--min-priority-bonus", "0.05",
             "--out-dir", v / "markers"])
        variants.append({"label": label, "score_threshold": float(score), "k": k, "nms_radius": nms,
                         "marker_path": str(v / "markers" / "markers_uint16.tif")})

    for score in ("0.25", "0.30", "0.35", "0.40", "0.45"):
        selector_variant(score, 3600, 6)
    for k in (900, 1800, 7200):
        selector_variant("0.35", k, 6)
    default_spatial = out / "score_0p35_k3600_n6" / "spatial" / "candidate_diagnostics_no_gt.csv"
    for nms in (4, 8, 10):
        markers(f"score_0p35_k3600_n{nms}", default_spatial, "0.35", 3600, nms)
    pd.DataFrame(variants).to_csv(out / "variants.csv", index=False)

    default_marker = out / "score_0p35_k3600_n6" / "markers" / "markers_uint16.tif"
    ws_args = ["--image", src["flat"], "--mask-marker", default_marker, "--mask-marker", src["baseline"],
               "--outside-policy", "drop", "--marker", f"baseline_existing={src['baseline']}",
               "--out-dir", out / "watershed"]
    for v in variants:
        ws_args += ["--marker", f"{v['label']}={v['marker_path']}"]
    run("watershed_fixed_mask", "e11_watershed", ws_args)
    results = evaluate(run, "evaluate_e09b", out / "watershed", src["gt"], "score_0p35_k3600_n6", out / "evaluation")
    merged = pd.DataFrame(variants).merge(results[results.evaluation_protocol == "seed_ignore"],
                                          left_on="label", right_on="marker_source", how="left")
    merged.to_csv(out / "development_sensitivity_results.csv", index=False)
    write_record(out, "phase1e_table4_sensitivity", run, t0, src, results, {"crop_id": DEV_CROP})
    print(merged[["label", "f1_iou050", "f1_iou075"]].to_string(index=False), flush=True)


def cropsize(root: Path) -> None:
    out = root / "crop_size"
    run, t0 = start(out)
    meta = json.loads((REPO_ROOT / "data/metadata/crops_sibgrapi.json").read_text(encoding="utf-8"))
    parents = ("EFRGP01_00_crop04", "EFRGP01_01_crop04")
    manifest = out / "e13_parent_manifest.csv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scan", "crop_id", "gray", "gt", "source_crop_dir"])
        writer.writeheader()
        for c in meta["crops"]:
            if c["crop_id"] in parents:
                writer.writerow({"scan": c["scan"], "crop_id": c["crop_id"], "gray": c["gray_crop_path"],
                                 "gt": c["gt_path"], "source_crop_dir": str(Path(c["gray_crop_path"]).parent)})
    run("e15_make_crop_size_suite", "e15_suite",
        ["--e13-manifest", manifest, "--e13-params-dir", RUNS / "_scale" / "auto_scale_params",
         "--suite-root", out / "suite", "--sizes", 128, 160, 200])
    subs = pd.read_csv(out / "suite" / "e15_crop_size_manifest.csv")
    for _, s in subs.iterrows():
        cmd = [sys.executable, "-I", "-m", "batos_v2.frozen", "--crop-id", s["crop_id"], "--runs-root", out / "e11",
               "--gray", s["gray"], "--gt", s["gt"],
               "--params-json", out / "suite" / "fixed_parent_params" / f"{s['crop_id']}_auto_params.json"]
        print(f"[cropsize] {s['crop_id']}", flush=True)
        done = subprocess.run([str(c) for c in cmd])
        if done.returncode != 0:
            raise RuntimeError(f"frozen run failed for {s['crop_id']}")
    rows = []
    for _, s in subs.iterrows():
        rec = json.loads((out / "e11" / s["crop_id"] / "frozen_run.json").read_text(encoding="utf-8"))
        for r in rec["e09b_results"]:
            rows.append({"crop_id": s["crop_id"], "parent_crop_id": s["parent_crop_id"],
                         "crop_size": s["crop_size"]} | r)
    results = pd.DataFrame(rows)
    results.to_csv(out / "crop_size_results.csv", index=False)
    write_record(out, "phase1e_table5_crop_size", run, t0, {"manifest": manifest}, results)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("what", choices=["mintree", "sensitivity", "cropsize"])
    parser.add_argument("--crop-id", nargs="*", default=[])
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args(argv)
    if args.what == "mintree":
        if not args.crop_id:
            parser.error("mintree needs --crop-id")
        for crop_id in args.crop_id:
            print(f"[mintree] {crop_id}", flush=True)
            mintree(crop_id, args.root)
    elif args.what == "sensitivity":
        sensitivity(args.root)
    else:
        cropsize(args.root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
