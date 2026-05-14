#!/usr/bin/env python3
import csv
from pathlib import Path
import statistics as st

inputs = [
    Path("results_batos/EFRGP01_00_multicrop/summary_multicrop.csv"),
    Path("results_batos/EFRGP01_00_extra_crops/summary_extra_crops.csv"),
]

rows = []

for path in inputs:
    if not path.exists():
        continue
    with open(path) as f:
        for r in csv.DictReader(f):
            r["gt_valid"] = int(r["gt_valid"])
            r["pred_labels"] = int(r["pred_labels"])
            for k in [
                "fg_iou",
                "fg_dice",
                "fg_precision",
                "fg_recall",
                "obj_recall_iou25",
                "obj_recall_iou50",
            ]:
                r[k] = float(r[k])
            rows.append(r)

selected = [r for r in rows if r["gt_valid"] >= 100]

out_dir = Path("results_batos/paper_tables")
out_dir.mkdir(parents=True, exist_ok=True)

with open(out_dir / "batos_selected_crops.csv", "w", newline="") as f:
    fieldnames = [
        "crop_id", "x", "y", "z", "size",
        "gt_valid", "pred_labels",
        "fg_iou", "fg_dice",
        "fg_precision", "fg_recall",
        "obj_recall_iou25", "obj_recall_iou50",
        "matched_iou50",
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in selected:
        w.writerow({k: r[k] for k in fieldnames})

summary_metrics = [
    "fg_iou",
    "fg_dice",
    "fg_precision",
    "fg_recall",
    "obj_recall_iou25",
    "obj_recall_iou50",
]

with open(out_dir / "batos_selected_crops_stats.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric", "mean", "median", "min", "max", "std"])
    for k in summary_metrics:
        vals = [r[k] for r in selected]
        std = st.stdev(vals) if len(vals) > 1 else 0.0
        w.writerow([
            k,
            f"{st.mean(vals):.6f}",
            f"{st.median(vals):.6f}",
            f"{min(vals):.6f}",
            f"{max(vals):.6f}",
            f"{std:.6f}",
        ])

print("selected crops:", len(selected))
print("wrote:", out_dir / "batos_selected_crops.csv")
print("wrote:", out_dir / "batos_selected_crops_stats.csv")
