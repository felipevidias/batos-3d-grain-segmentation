#!/usr/bin/env python3
"""Phase 1 closure: aggregate the 12 frozen batos-paper-v2 runs into results/batos_sibgrapi/.

Per crop (from C:\\BA_TOS_DATA\\runs\\batos_paper_v2_frozen\\<crop_id>):
  1. e09b (frozen, unmodified) + IoU 0.25 via scripts/e09b_with_iou025.py, for
     BA-TOS_ToS (e11/watershed: tos_e11_min), the Otsu+WS control (edt_otsu) and
     BA-TOS_sel (e23c/watershed: batos_selective).
  2. Same evaluation after the legacy prediction size filter (recovered
     filter_predicted_labels_by_size.py --min-size 500, verbatim, IDs preserved), labelled size_filter=500.
  3. Check that 0.50/0.75 from the wrapper equal the plain e09b run, and that ToS and sel segmentations are
     bitwise identical on gate-negative crops.
Outputs: per_crop_metrics.csv, crop_macro_summary.csv, pooled_object_summary.csv, table2_comparison.csv,
REPORT.md, and results/all_runs.csv (plan columns, runtime and peak memory).
Primary protocol = seed_ignore (paper Sec. III-C). Pooled and crop-macro tables are kept separate.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

REPO = Path(__file__).resolve().parents[1]
RUNS = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen")
RERUN = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen_rerun")  # determinism reruns (4 crops)
MEMRUN = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen_mem")  # full reruns with process-tree peak memory
OUT = REPO / "results" / "batos_sibgrapi"
E09B = Path(r"C:\BA_TOS_DATA\linux_import\chat_recovered\e13_crop05_metadata_full_20260725_182408\audit"
            r"\e09b_ignore_aware_audit.py")
FILTER = Path(r"C:\BA_TOS_DATA\linux_import\zsh_recovered"
              r"\code__segmentation__batos__filter_predicted_labels_by_size.py")
WRAPPER = REPO / "scripts" / "e09b_with_iou025.py"
CROPS = [f"EFRGP01_0{s}_crop0{c}" for s in (0, 1) for c in (2, 3, 4, 5, 6, 7)]
ROUTES = {  # route label -> (watershed dir relative to crop, configuration, marker label)
    "BA-TOS_ToS": ("e11/watershed", "tos_e11_min_distance_c1", "tos_e11_min"),
    "Otsu+WS(edt_otsu)": ("e11/watershed", "edt_otsu_distance_c1", "edt_otsu"),
    "BA-TOS_sel": ("e23c/watershed", "batos_selective_distance_c1", "batos_selective"),
}
PAPER_TABLE2 = {  # crop-macro, Table II of the camera-ready paper
    "BA-TOS_ToS": {"mean": 0.824, "sd": 0.367, "median": 0.985, "min": 0.012, "mean_075": 0.817},
    "BA-TOS_sel": {"mean": 0.906, "sd": 0.210, "median": 0.985, "min": 0.260, "mean_075": 0.899},
    "Otsu+WS(edt_otsu)": {"mean": 0.830, "sd": 0.182, "median": 0.882, "min": 0.260, "mean_075": 0.770},
}
PY = sys.executable


def call(args: list) -> None:
    done = subprocess.run([PY, "-I", *[str(a) for a in args]], capture_output=True, text=True)
    if done.returncode != 0:
        raise RuntimeError(f"failed: {args[:2]}\n{done.stderr[-1500:]}")


def array_sha(path: Path) -> str:
    return hashlib.sha256(np.ascontiguousarray(tifffile.imread(path)).tobytes()).hexdigest()


def evaluate(crop: str, route: str, size_filter: str) -> pd.DataFrame:
    crop_dir = RUNS / crop
    ws_rel, configuration, label = ROUTES[route]
    ws_dir = crop_dir / ws_rel
    rows = pd.read_csv(ws_dir / "watershed_results.csv")
    row = rows[rows["configuration"] == configuration]
    if len(row) != 1:
        raise RuntimeError(f"{crop}: configuration {configuration} not found in {ws_dir}")
    tag = f"{label}_{size_filter}"
    work = crop_dir / "phase1_closure" / tag
    if work.exists():
        shutil.rmtree(work)
    e09_dir = work / "e09_dir"
    (e09_dir / "segmentations").mkdir(parents=True)
    seg = ws_dir / "segmentations" / f"{configuration}.tif"
    if size_filter == "none":
        shutil.copy2(seg, e09_dir / "segmentations" / f"{configuration}.tif")
    else:
        # No --relabel: e09b's seed/overlap-ignore treats prediction label == marker label (the watershed
        # keeps marker IDs). Sequential relabelling breaks that link and corrupts the ignore sets (observed
        # 2026-11-04: ToS crop04 0.985 -> 0.726 with only 4 objects removed). Without it the filter only
        # zeroes objects < min-size and keeps all other IDs.
        call([FILTER, "--labels", seg, "--out-dir", work / "filter", "--prefix", configuration,
              "--min-size", size_filter])
        shutil.copy2(work / "filter" / f"{configuration}_min{size_filter}_labels.tif",
                     e09_dir / "segmentations" / f"{configuration}.tif")
    row.to_csv(e09_dir / "watershed_results.csv", index=False)
    marker = row["normalized_marker"].iloc[0]
    call([WRAPPER, E09B, "--e09-dir", e09_dir, "--gt", Path(rf"C:\BA_TOS_DATA\crops\{crop}\labels.tif"),
          "--marker", f"{label}={marker}", "--tos-label", label, "--out-dir", work / "evaluation"])
    res = pd.read_csv(work / "evaluation" / "ignore_aware_results.csv")
    res.insert(0, "size_filter", size_filter)
    res.insert(0, "route", route)
    res.insert(0, "crop_id", crop)
    return res


def run_record(crop: str) -> dict:
    """Runtime + provenance from the primary run (sequential, uncontended); peak memory only from the
    memory-run root, whose records carry peak_memory_method (process-tree peak). Earlier records measured
    only the ~5 MB venv launcher and are never used for memory."""
    rec = json.loads((RUNS / crop / "frozen_run.json").read_text(encoding="utf-8"))
    rec = rec | {"_source": str(RUNS / crop / "frozen_run.json"), "peak_memory_mb": None, "_mem_source": None}
    mem_path = MEMRUN / crop / "frozen_run.json"
    if mem_path.exists():
        mem = json.loads(mem_path.read_text(encoding="utf-8"))
        if all("peak_memory_method" in s for s in mem.get("steps", [])) and mem.get("steps"):
            rec["peak_memory_mb"] = mem["peak_memory_mb"]
            rec["_mem_source"] = str(mem_path)
    return rec


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    frames, checks = [], []
    for crop in CROPS:
        if not (RUNS / crop / "frozen_run.json").exists():
            raise FileNotFoundError(f"missing frozen run for {crop}")
        for route in ROUTES:
            for size_filter in ("none", "500"):
                frames.append(evaluate(crop, route, size_filter))
                print(f"[agg] {crop} {route} filter={size_filter}", flush=True)
        # check 1: wrapper 0.50/0.75 == plain e09b run
        plain = pd.read_csv(RUNS / crop / "e11" / "evaluation" / "ignore_aware_results.csv")
        mine = frames[-6]  # BA-TOS_ToS, none
        a = plain[(plain.marker_source == "tos_e11_min")].set_index("evaluation_protocol")
        b = mine.set_index("evaluation_protocol")
        cols = ["tp_iou050", "fp_iou050", "fn_iou050", "f1_iou050", "tp_iou075", "f1_iou075"]
        same_metrics = bool((a.loc[b.index, cols].to_numpy() == b[cols].to_numpy()).all())
        # check 2: ToS vs sel bitwise when the gate is negative
        route = json.loads((RUNS / crop / "e23c" / "route" / "route_manifest.json").read_text(encoding="utf-8"))
        gate = route.get("gate", route)
        tos_sha = array_sha(RUNS / crop / "e11/watershed/segmentations/tos_e11_min_distance_c1.tif")
        sel_sha = array_sha(RUNS / crop / "e23c/watershed/segmentations/batos_selective_distance_c1.tif")
        checks.append({"crop_id": crop, "wrapper_050_075_equal_plain_e09b": same_metrics,
                       "gate_triggered": bool(gate["gate_triggered"]),
                       "region_voxels_median": gate["region_voxels_median"],
                       "region_voxels_q75": gate["region_voxels_q75"],
                       "tos_sel_segmentation_bitwise_identical": tos_sha == sel_sha})
    per = pd.concat(frames, ignore_index=True)
    per["overseg_ratio"] = per["prediction_count_evaluable"] / per["valid_gt_count"]
    per.to_csv(OUT / "per_crop_metrics_all_protocols.csv", index=False)
    seed = per[per.evaluation_protocol == "seed_ignore"].copy()
    keep = ["crop_id", "route", "size_filter", "valid_gt_count", "prediction_count_evaluable", "overseg_ratio"]
    for s in ("025", "050", "075"):
        keep += [f"tp_iou{s}", f"fp_iou{s}", f"fn_iou{s}", f"precision_iou{s}", f"recall_iou{s}", f"f1_iou{s}"]
    seed[keep].to_csv(OUT / "per_crop_metrics_seed_ignore.csv", index=False)
    pd.DataFrame(checks).to_csv(OUT / "checks.csv", index=False)

    macro_rows, pooled_rows = [], []
    for (route, size_filter), g in seed.groupby(["route", "size_filter"], sort=False):
        m = {"summary_type": "crop_macro", "route": route, "size_filter": size_filter, "n_crops": len(g)}
        p = {"summary_type": "pooled_object", "route": route, "size_filter": size_filter, "n_crops": len(g),
             "valid_gt": int(g.valid_gt_count.sum()), "pred_eval": int(g.prediction_count_evaluable.sum())}
        p["overseg_ratio"] = p["pred_eval"] / p["valid_gt"]
        for s in ("025", "050", "075"):
            f1 = g[f"f1_iou{s}"].astype(float)
            m.update({f"f1_{s}_mean": f1.mean(), f"f1_{s}_sd_sample": f1.std(ddof=1),
                      f"f1_{s}_sd_population": f1.std(ddof=0), f"f1_{s}_median": f1.median(),
                      f"f1_{s}_min": f1.min()})
            tp, fp, fn = (int(g[f"{k}_iou{s}"].sum()) for k in ("tp", "fp", "fn"))
            pr = tp / (tp + fp) if tp + fp else 0.0
            rc = tp / (tp + fn) if tp + fn else 0.0
            p.update({f"tp_{s}": tp, f"fp_{s}": fp, f"fn_{s}": fn, f"precision_{s}": pr, f"recall_{s}": rc,
                      f"f1_{s}": 2 * pr * rc / (pr + rc) if pr + rc else 0.0})
        m["overseg_ratio_mean"] = float(g.overseg_ratio.mean())
        macro_rows.append(m)
        pooled_rows.append(p)
    macro, pooled = pd.DataFrame(macro_rows), pd.DataFrame(pooled_rows)
    macro.to_csv(OUT / "crop_macro_summary.csv", index=False)
    pooled.to_csv(OUT / "pooled_object_summary.csv", index=False)

    comp = []
    for route, paper in PAPER_TABLE2.items():
        row = macro[(macro.route == route) & (macro.size_filter == "none")].iloc[0]
        for key, ours in (("mean", row.f1_050_mean), ("sd", row.f1_050_sd_sample), ("median", row.f1_050_median),
                          ("min", row.f1_050_min), ("mean_075", row.f1_075_mean)):
            comp.append({"route": route, "statistic": key, "paper": paper[key], "ours": round(float(ours), 4),
                         "diff": round(float(ours) - paper[key], 4)})
        comp.append({"route": route, "statistic": "sd_population", "paper": paper["sd"],
                     "ours": round(float(row.f1_050_sd_population), 4),
                     "diff": round(float(row.f1_050_sd_population) - paper["sd"], 4)})
    comparison = pd.DataFrame(comp)
    comparison.to_csv(OUT / "table2_comparison.csv", index=False)

    # all_runs.csv (plan columns); previous file (Phase 1c port runs) is kept as all_runs_phase1c_port.csv
    all_runs = REPO / "results" / "all_runs.csv"
    if all_runs.exists() and not (REPO / "results" / "all_runs_phase1c_port.csv").exists():
        all_runs.rename(REPO / "results" / "all_runs_phase1c_port.csv")
    rows = []
    for _, r in seed.iterrows():
        rec = run_record(r.crop_id)
        rows.append({
            "experiment_id": f"batos_paper_v2_frozen/{r.crop_id}", "method": r.route,
            "model": "batos-paper-v2-frozen (E-series verbatim)", "volume": r.crop_id[:10], "crop": r.crop_id,
            "device": "cpu (AMD Ryzen 7 5700X)",
            "parameters": f"size_filter={r.size_filter};protocol=seed_ignore;K=3600;tau_s=0.35;NMS=6",
            "runtime_seconds": rec.get("seconds_total"), "peak_memory_mb": rec.get("peak_memory_mb"),
            "tp": r.tp_iou050, "fp": r.fp_iou050, "fn": r.fn_iou050, "precision": r.precision_iou050,
            "recall": r.recall_iou050, "f1_025": r.f1_iou025, "f1_050": r.f1_iou050, "f1_075": r.f1_iou075,
            "ap50": None, "ap75": None, "mean_matched_iou": r.get("mean_matched_iou050"),
            "n_pred": r.prediction_count_evaluable, "n_gt": r.valid_gt_count,
            "git_commit": rec["git"]["commit"], "run_json": rec["_source"], "memory_run_json": rec["_mem_source"],
        })
    pd.DataFrame(rows).to_csv(all_runs, index=False)
    print(comparison.to_string(index=False))
    print(pd.DataFrame(checks).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
