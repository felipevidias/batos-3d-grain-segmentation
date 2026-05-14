#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT = Path(".")
IN = ROOT / "results_batos/paper_tables/grain_centered_analysis/grain_centered_detail_per_gt.csv"
OUT = ROOT / "results_batos/paper_tables/grain_centered_analysis"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(IN)

# Outliers do BA-TOS
batos = df[df["method"] == "BA-TOS"].copy()

# Casos ruins por IoU
low_iou = batos.sort_values("best_iou", ascending=True).head(50)
low_iou.to_csv(OUT / "batos_worst_50_by_iou.csv", index=False)

# Casos ruins por erro de centroide
centroid_out = batos.sort_values("centroid_error", ascending=False).head(50)
centroid_out.to_csv(OUT / "batos_worst_50_by_centroid_error.csv", index=False)

# Casos ruins por erro de volume
vol_out = batos.sort_values("relative_volume_error", ascending=False).head(50)
vol_out.to_csv(OUT / "batos_worst_50_by_volume_error.csv", index=False)

# Comparação por crop: onde BA-TOS sofre mais?
crop_summary = (
    batos.groupby(["scan", "crop_id"], as_index=False)
    .agg(
        n_grains=("gt_label", "count"),
        rec050=("recovered_iou050", "mean"),
        rec075=("recovered_iou075", "mean"),
        median_iou=("best_iou", "median"),
        mean_iou=("best_iou", "mean"),
        median_centroid_error=("centroid_error", "median"),
        mean_centroid_error=("centroid_error", "mean"),
        median_volume_error=("relative_volume_error", "median"),
        mean_volume_error=("relative_volume_error", "mean"),
    )
    .sort_values(["rec075", "median_iou"], ascending=[True, True])
)

crop_summary.to_csv(OUT / "batos_grain_centered_summary_by_crop.csv", index=False)

print("wrote:", OUT / "batos_worst_50_by_iou.csv")
print("wrote:", OUT / "batos_worst_50_by_centroid_error.csv")
print("wrote:", OUT / "batos_worst_50_by_volume_error.csv")
print("wrote:", OUT / "batos_grain_centered_summary_by_crop.csv")

print()
print("=== BA-TOS crops sorted by weakest Rec@.75 / median IoU ===")
print(crop_summary.to_string(index=False))
