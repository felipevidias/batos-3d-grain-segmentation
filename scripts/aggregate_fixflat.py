#!/usr/bin/env python3
"""Aggregate batos-paper-v2-fixflat (12 crops) into results/batos_v2_fixflat/.

batos-paper-v2-fixflat = the frozen batos-paper-v2 pipeline with one change: the lambda=1 selector flatten takes its
offset over the specimen (bg_ref = median(blur >= 1), batos_v2.flatten_fixed --variant A) instead of
median(blur > 0). Scale and openrec are read from the frozen run; every later stage is the frozen call.
Runs: C:\\BA_TOS_DATA\\runs\\batos_paper_v2_fixflat\\A\\<crop> (driver: python -I -m batos_v2.flatten_fix_run).

Evaluation is the same as the Phase 1 closure (scripts/aggregate_sibgrapi.py, whose evaluate() is reused): frozen
e09b + IoU 0.25 wrapper, seed-ignore primary protocol, size filter none and 500 (no relabel). Outputs:
per_crop_metrics_seed_ignore.csv, per_crop_metrics_all_protocols.csv, crop_macro_summary.csv,
pooled_object_summary.csv, table2_fixflat_vs_frozen.csv, gate_and_runtime.csv, checks.csv; rows upserted into
results/all_runs.csv (experiment_id batos_paper_v2_fixflat/<crop>, key experiment_id+method+parameters).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.extinction_run import upsert  # noqa: E402

spec = importlib.util.spec_from_file_location("agg_frozen", REPO / "scripts" / "aggregate_sibgrapi.py")
agg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(agg)

RUNS = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_fixflat\A")
FROZEN = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen")
OUT = REPO / "results" / "batos_v2_fixflat"
FROZEN_TABLES = REPO / "results" / "batos_sibgrapi"
ROUTES = ["BA-TOS_ToS", "BA-TOS_sel"]  # the driver runs no Otsu+WS control (count-matched EDT-Otsu is only the gate fallback)


def summaries(seed: pd.DataFrame):
    macro_rows, pooled_rows = [], []
    for (route, size_filter), g in seed.groupby(["route", "size_filter"], sort=False):
        m = {"summary_type": "crop_macro", "route": route, "size_filter": size_filter, "n_crops": len(g)}
        p = {"summary_type": "pooled_object", "route": route, "size_filter": size_filter, "n_crops": len(g),
             "valid_gt": int(g.valid_gt_count.sum()), "pred_eval": int(g.prediction_count_evaluable.sum())}
        for s in ("025", "050", "075"):
            f1 = g[f"f1_iou{s}"].astype(float)
            m.update({f"f1_{s}_mean": f1.mean(), f"f1_{s}_sd_sample": f1.std(ddof=1), f"f1_{s}_median": f1.median(),
                      f"f1_{s}_min": f1.min(), f"f1_{s}_argmin": g.loc[f1.idxmin(), "crop_id"]})
            tp, fp, fn = (int(g[f"{k}_iou{s}"].sum()) for k in ("tp", "fp", "fn"))
            pr = tp / (tp + fp) if tp + fp else 0.0
            rc = tp / (tp + fn) if tp + fn else 0.0
            p.update({f"tp_{s}": tp, f"fp_{s}": fp, f"fn_{s}": fn, f"precision_{s}": pr, f"recall_{s}": rc,
                      f"f1_{s}": 2 * pr * rc / (pr + rc) if pr + rc else 0.0})
        macro_rows.append(m)
        pooled_rows.append(p)
    return pd.DataFrame(macro_rows), pd.DataFrame(pooled_rows)


def main() -> int:
    agg.RUNS = RUNS  # evaluate() reads <RUNS>/<crop>/{e11,e23c}/watershed and writes <RUNS>/<crop>/phase1_closure
    OUT.mkdir(parents=True, exist_ok=True)
    frames, checks, gates = [], [], []
    for crop in agg.CROPS:
        if not (RUNS / crop / "DONE").exists():
            raise FileNotFoundError(f"run not complete: {RUNS / crop}")
        for route in ROUTES:
            for size_filter in ("none", "500"):
                frames.append(agg.evaluate(crop, route, size_filter))
        print(f"[agg] {crop}", flush=True)
        rec = json.loads((RUNS / crop / "run.json").read_text(encoding="utf-8"))
        frec = json.loads((FROZEN / crop / "frozen_run.json").read_text(encoding="utf-8"))
        upstream = sum(s["seconds"] for s in frec["steps"] if s.get("step") in ("scale", "openrec_v2"))
        fixed = sum(s["seconds"] for s in rec["steps"])
        route = json.loads((RUNS / crop / "e23c" / "route" / "route_manifest.json").read_text(encoding="utf-8"))
        gate = route.get("gate", route)
        fgate = json.loads((FROZEN / crop / "e23c" / "route" / "route_manifest.json").read_text(encoding="utf-8"))
        fgate = fgate.get("gate", fgate)
        sha = rec["segmentation_sha256"]
        gates.append({"crop_id": crop, "bg_ref_fixflat": rec["bg_ref"], "gate_triggered": bool(gate["gate_triggered"]),
                      "region_voxels_median": gate["region_voxels_median"],
                      "region_voxels_q75": gate["region_voxels_q75"],
                      "frozen_gate_triggered": bool(fgate["gate_triggered"]),
                      "frozen_region_voxels_median": fgate["region_voxels_median"],
                      "frozen_region_voxels_q75": fgate["region_voxels_q75"],
                      "tos_sel_bitwise_identical": sha["BA-TOS_ToS"] == sha["BA-TOS_sel"],
                      "seconds_flatten_onward": round(fixed, 1), "seconds_upstream_frozen_scale_openrec": upstream,
                      "seconds_end_to_end": round(fixed + upstream, 1),
                      "frozen_seconds_total": frec.get("seconds_total"),
                      "peak_stage_memory_mb": max(s.get("peak_memory_mb", 0) for s in rec["steps"])})
    per = pd.concat(frames, ignore_index=True)
    per["overseg_ratio"] = per["prediction_count_evaluable"] / per["valid_gt_count"]
    per.to_csv(OUT / "per_crop_metrics_all_protocols.csv", index=False)
    seed = per[per.evaluation_protocol == "seed_ignore"].copy()
    keep = ["crop_id", "route", "size_filter", "valid_gt_count", "prediction_count_evaluable", "overseg_ratio"]
    for s in ("025", "050", "075"):
        keep += [f"tp_iou{s}", f"fp_iou{s}", f"fn_iou{s}", f"precision_iou{s}", f"recall_iou{s}", f"f1_iou{s}"]
    seed[keep].to_csv(OUT / "per_crop_metrics_seed_ignore.csv", index=False)
    gate_df = pd.DataFrame(gates)
    gate_df.to_csv(OUT / "gate_and_runtime.csv", index=False)

    # check: e09b (this evaluation) vs the driver's evaluate_v2 at IoU 0.50 / 0.75 (no filter)
    for crop in agg.CROPS:
        v2 = pd.read_csv(RUNS / crop / "metrics.csv").set_index("route")
        for route in ("BA-TOS_ToS", "BA-TOS_sel"):
            a = seed[(seed.crop_id == crop) & (seed.route == route) & (seed.size_filter == "none")].iloc[0]
            b = v2.loc[route]
            checks.append({"crop_id": crop, "route": route,
                           "e09b_equals_evaluate_v2": all(int(a[f"{k}_iou{s}"]) == int(b[f"{k}_iou{s}"])
                                                          for k in ("tp", "fp", "fn") for s in ("050", "075"))})
    pd.DataFrame(checks).to_csv(OUT / "checks.csv", index=False)

    macro, pooled = summaries(seed)
    macro.to_csv(OUT / "crop_macro_summary.csv", index=False)
    pooled.to_csv(OUT / "pooled_object_summary.csv", index=False)
    fmacro = pd.read_csv(FROZEN_TABLES / "crop_macro_summary.csv")
    comp = []
    for route in ROUTES:
        for sf in ("none", "500"):
            a = macro[(macro.route == route) & (macro.size_filter.astype(str) == sf)].iloc[0]
            b = fmacro[(fmacro.route == route) & (fmacro.size_filter.astype(str) == sf)].iloc[0]
            for st in ("f1_050_mean", "f1_050_sd_sample", "f1_050_median", "f1_050_min", "f1_075_mean"):
                comp.append({"route": route, "size_filter": sf, "statistic": st, "frozen": round(float(b[st]), 4),
                             "fixflat": round(float(a[st]), 4), "diff": round(float(a[st]) - float(b[st]), 4)})
    pd.DataFrame(comp).to_csv(OUT / "table2_fixflat_vs_frozen.csv", index=False)

    git = json.loads((RUNS / agg.CROPS[0] / "run.json").read_text(encoding="utf-8"))["git"]["commit"]
    rows = []
    g = gate_df.set_index("crop_id")
    for _, r in seed.iterrows():
        rows.append({
            "experiment_id": f"batos_paper_v2_fixflat/{r.crop_id}", "method": r.route,
            "model": "batos-paper-v2-fixflat (frozen + flatten bg_ref = median(blur>=1))",
            "volume": r.crop_id[:10], "crop": r.crop_id, "device": "cpu (AMD Ryzen 7 5700X)",
            "parameters": f"size_filter={r.size_filter};protocol=seed_ignore;K=3600;tau_s=0.35;NMS=6;flatten_bg_ref=A",
            "runtime_seconds": g.loc[r.crop_id, "seconds_end_to_end"],
            "peak_memory_mb": g.loc[r.crop_id, "peak_stage_memory_mb"],
            "tp": r.tp_iou050, "fp": r.fp_iou050, "fn": r.fn_iou050, "precision": r.precision_iou050,
            "recall": r.recall_iou050, "f1_025": r.f1_iou025, "f1_050": r.f1_iou050, "f1_075": r.f1_iou075,
            "ap50": None, "ap75": None, "n_pred": r.prediction_count_evaluable, "n_gt": r.valid_gt_count,
            "git_commit": git, "run_json": str(RUNS / r.crop_id / "run.json")})
    upsert(REPO / "results" / "all_runs.csv", rows, ["experiment_id", "method", "parameters"])
    print(pd.DataFrame(comp).to_string(index=False))
    print(gate_df.to_string(index=False))
    print(pd.DataFrame(checks).e09b_equals_evaluate_v2.value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
