#!/usr/bin/env python3
"""Phase 2: cpsam_v2 on the synthetic scale set (factors 0.7 / 1.4 / 2.0 of the 12 dev crops), vs frozen BA-TOS v3.

Synthetic crops: scripts/make_synthetic_scale.py (C:\\BA_TOS_DATA\\crops_synth\\x<f>\\<crop>\\{gray,labels}.tif,
manifest data/metadata/crops_manifest_synth.csv). Valid GT threshold min_gt_voxels = round(1000 f^3) from that manifest
(same physical grains), face-touching GT excluded, as in the v3 scale test.
Cellpose predictions: C:\\BA_TOS_DATA\\runs\\cellpose\\synth_cpsam_v2\\<crop>\\labels.tif (run_cellpose.py, fixed
cpsam_v2 config, no diameter, no tuning). Evaluation: evaluate_v2 no-marker mode (fraction 0.50), size filter none
(primary) and 500 f^3 (secondary, labelled; uses the known synthetic factor, the analogue of v3's own 500 s^3 filter).
v3 rows: read from the v3 scale test (C:\\BA_TOS_DATA\\runs\\v3_scale_test\\eval, method v3, route BA-TOS_sel,
mode no_markers_fraction050; size filter none / own). Factor 1.0 = the original dev crops (cpsam_v2 eval cache).
Batos env. Writes results/cellpose/synthetic_scale/{per_crop.csv, crop_macro.csv, pooled.csv, tables.md} and upserts
all_runs rows cellpose_synth_cpsam_v2/<crop>.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from batos_v2.evaluate_v2 import evaluate  # noqa: E402
from evaluate_cellpose import EVAL_KW, MODE, RUNS, V3_DEV_EVAL, git_state, size_filter, upsert  # noqa: E402

TAG = "synth_cpsam_v2"
SYNTH = Path(r"C:\BA_TOS_DATA\crops_synth")
OUT = REPO / "results" / "cellpose" / "synthetic_scale"
TS = ("025", "050", "075")


def fdir(f: float) -> str:
    return "x" + (f"{f:g}".replace(".", "p"))


def cellpose_rows(c) -> pd.DataFrame:
    cache = RUNS / TAG / "eval" / f"{c.crop_id}.csv"
    if cache.exists():
        return pd.read_csv(cache, dtype={"size_filter": str})
    d = RUNS / TAG / c.crop_id
    gt = tifffile.imread(SYNTH / fdir(c.factor) / c.crop_id / "labels.tif")
    if (d / "FAILED.json").exists():
        pred, status = np.zeros(gt.shape, np.uint16), "failed"
    elif (d / "labels.tif").exists() and (d / "run.json").exists():
        pred, status = tifffile.imread(d / "labels.tif"), "ok"
    else:
        raise SystemExit(f"missing prediction {d}")
    if pred.shape != gt.shape:
        raise SystemExit(f"shape mismatch {c.crop_id}: {pred.shape} vs {gt.shape}")
    rows = []
    for flt, ms in (("none", None), ("own", round(500 * c.factor ** 3))):
        m, _, _ = evaluate(pred if ms is None else size_filter(pred, ms), gt, None,
                           min_gt_voxels=int(c.min_gt_voxels), **EVAL_KW)
        rows.append({"crop_id": c.crop_id, "dev_crop_id": c.dev_crop_id, "factor": c.factor,
                     "method": "cellpose:cpsam_v2", "status": status, "size_filter": flt, "size_filter_voxels": ms,
                     "mode": MODE, **m})
    df = pd.DataFrame(rows)
    cache.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False)
    return df


def main() -> int:
    man = pd.read_csv(REPO / "data" / "metadata" / "crops_manifest_synth.csv")
    parts = [cellpose_rows(c) for c in man.itertuples()]
    dev = sorted(man.dev_crop_id.unique())
    for crop in dev:  # factor 1.0 = original dev crops
        d = pd.read_csv(RUNS / "cpsam_v2" / "eval" / f"{crop}.csv", dtype={"size_filter": str})
        d = d.assign(dev_crop_id=crop, factor=1.0, size_filter=d.size_filter.replace({"500": "own"}))
        parts.append(d)
    v3 = []
    for p in V3_DEV_EVAL.glob("*.csv"):
        d = pd.read_csv(p, dtype={"size_filter": str})
        d = d[(d.method == "v3") & (d.route == "BA-TOS_sel") & (d["mode"] == MODE)]
        v3.append(d.assign(method="v3", status="ok"))
    v3 = pd.concat(v3, ignore_index=True)
    v3 = v3[v3.dev_crop_id.isin(dev)]
    per = pd.concat(parts + [v3], ignore_index=True)
    keep = ["crop_id", "dev_crop_id", "factor", "method", "status", "size_filter", "mode", "valid_gt_count",
            "prediction_count_evaluable", "mean_matched_iou050"] + [f"{k}_iou{s}" for s in TS
                                                                    for k in ("tp", "fp", "fn", "f1")]
    per = per[keep].sort_values(["size_filter", "factor", "method", "crop_id"])
    for (flt, f, m), g in per.groupby(["size_filter", "factor", "method"]):
        if len(g) != 12:
            raise SystemExit(f"expected 12 crops for {flt} {f} {m}, got {len(g)}")
    macro, pooled = [], []
    for (flt, f, m), g in per.groupby(["size_filter", "factor", "method"]):
        b = {"size_filter": flt, "factor": f, "method": m, "n": len(g)}
        mm = dict(b)
        for s in TS:
            mm |= {f"f1_{s}_mean": g[f"f1_iou{s}"].mean(), f"f1_{s}_sd": g[f"f1_iou{s}"].std(ddof=1),
                   f"f1_{s}_min": g[f"f1_iou{s}"].min()}
        mm |= {"argmin_050": g.loc[g.f1_iou050.idxmin(), "crop_id"], "miou_050": g.mean_matched_iou050.mean(),
               "n_failed": int((g.status != "ok").sum())}
        macro.append(mm)
        pp = dict(b) | {"gt": int(g.valid_gt_count.sum()), "pred": int(g.prediction_count_evaluable.sum())}
        for s in TS:
            tp, fp, fn = (int(g[f"{k}_iou{s}"].sum()) for k in ("tp", "fp", "fn"))
            pp |= {f"tp_{s}": tp, f"fp_{s}": fp, f"fn_{s}": fn, f"f1_{s}": 2 * tp / (2 * tp + fp + fn)}
        pooled.append(pp)
    macro, pooled = pd.DataFrame(macro), pd.DataFrame(pooled)
    OUT.mkdir(parents=True, exist_ok=True)
    per.to_csv(OUT / "per_crop.csv", index=False)
    macro.to_csv(OUT / "crop_macro.csv", index=False)
    pooled.to_csv(OUT / "pooled.csv", index=False)

    L = ["# Synthetic scale set: cpsam_v2 vs BA-TOS v3 (sel), no-marker mode (fraction 0.50)", "",
         "12 dev crops resampled by f about their centre; valid GT >= round(1000 f^3) voxels, not face-touching. "
         "f = 1 is the original dev crops. Size filter `own` = 500 f^3 for cpsam_v2, 500 s^3 (v3's GT-free scale) for v3.",
         ""]
    for flt in ("none", "own"):
        L += [f"## Crop-macro, size filter {flt}", "",
              "| factor | method | n | F1@0.50 mean ± SD | min (crop) | F1@0.25 | F1@0.75 | mIoU@0.50 | failed |",
              "|---|---|---|---|---|---|---|---|---|"]
        for r in macro[macro.size_filter == flt].itertuples():
            L.append(f"| {r.factor:g} | {r.method} | {r.n} | {r.f1_050_mean:.3f} ± {r.f1_050_sd:.3f} | "
                     f"{r.f1_050_min:.3f} ({r.argmin_050}) | {r.f1_025_mean:.3f} | {r.f1_075_mean:.3f} | "
                     f"{r.miou_050:.3f} | {r.n_failed} |")
        L += ["", f"## Pooled, size filter {flt}", "",
              "| factor | method | GT | pred | TP@0.50 | FP | FN | F1@0.25 | F1@0.50 | F1@0.75 |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for r in pooled[pooled.size_filter == flt].itertuples():
            L.append(f"| {r.factor:g} | {r.method} | {r.gt} | {r.pred} | {r.tp_050} | {r.fp_050} | {r.fn_050} | "
                     f"{r.f1_025:.3f} | {r.f1_050:.3f} | {r.f1_075:.3f} |")
        L.append("")
    (OUT / "tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    git = git_state()["commit"]
    rows = []
    for r in per[(per.method == "cellpose:cpsam_v2") & (per.factor != 1.0)].itertuples():
        rj = RUNS / TAG / r.crop_id / ("run.json" if r.status == "ok" else "FAILED.json")
        rec = json.loads(rj.read_text(encoding="utf-8"))
        rows.append({"experiment_id": f"cellpose_{TAG}/{r.crop_id}", "method": "cellpose", "model": rec["model"],
                     "volume": r.dev_crop_id[:10], "crop": r.crop_id, "device": rec["environment"]["device"],
                     "parameters": f"size_filter={r.size_filter};eval=evaluate_v2:{MODE};split=synthetic_scale;"
                                   f"factor={r.factor:g};tag={TAG}",
                     "runtime_seconds": rec.get("runtime_seconds"), "peak_memory_mb": rec.get("peak_host_rss_mb"),
                     "tp": r.tp_iou050, "fp": r.fp_iou050, "fn": r.fn_iou050, "f1_025": r.f1_iou025,
                     "f1_050": r.f1_iou050, "f1_075": r.f1_iou075, "mean_matched_iou": r.mean_matched_iou050,
                     "n_pred": r.prediction_count_evaluable, "n_gt": r.valid_gt_count, "git_commit": git,
                     "run_json": str(rj)})
    upsert(REPO / "results" / "all_runs.csv", rows, ["experiment_id", "method", "parameters"])
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
