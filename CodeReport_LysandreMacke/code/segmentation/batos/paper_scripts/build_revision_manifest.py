#!/usr/bin/env python3
from pathlib import Path
import csv
import re

ROOT = Path(".")
OUT = Path("results_batos/paper_tables/final_revision")
OUT.mkdir(parents=True, exist_ok=True)

# crops usados no artigo
CROPS = {
    "crop02": ("600", "800", "300"),
    "crop03": ("900", "600", "300"),
    "crop04": ("500", "500", "700"),
    "crop05": ("800", "900", "700"),
    "crop06": ("350", "650", "500"),
    "crop07": ("700", "400", "500"),
}

SCANS = ["EFRGP01_00", "EFRGP01_01"]

def find_one(patterns):
    hits = []
    for pat in patterns:
        hits.extend(ROOT.glob(pat))
    hits = [p for p in hits if p.exists() and p.is_file()]
    # Prefer min500/final when available
    hits = sorted(
        set(hits),
        key=lambda p: (
            "min500" not in str(p),
            "pipeline_final" not in str(p),
            len(str(p))
        )
    )
    return hits[0] if hits else None

def scan_base(scan):
    if scan == "EFRGP01_00":
        return "results_batos/EFRGP01_00_multicrop"
    return "results_batos/EFRGP01_01_selected_crops"

rows = []

for scan in SCANS:
    base = scan_base(scan)

    for crop_id, (x, y, z) in CROPS.items():
        # GT label
        gt_patterns = [
            f"{base}/{crop_id}_x{x}_y{y}_z{z}_s200/00_input/*label_crop.tif",
            f"{base}/{crop_id}_x{x}_y{y}_z{z}_s200/00_input/*_label_crop.tif",
            f"results_batos/{scan}_test_{crop_id}/00_input/*label_crop.tif",
            f"results_batos/**/{scan}*{crop_id}*/00_input/*label_crop.tif",
            f"results_batos/**/{crop_id}_x{x}_y{y}_z{z}_s200_label_crop.tif",
        ]
        gt = find_one(gt_patterns)

        # BA-TOS final min500
        batos_patterns = [
            f"{base}/{crop_id}_x{x}_y{y}_z{z}_s200/batos_pipeline_final_min500/05_watershed/*min500*labels.tif",
            f"{base}/{crop_id}_x{x}_y{y}_z{z}_s200/batos_pipeline_final_min500/**/*.tif",
            f"results_batos/paper_tables/batos_size_filter_sweep/min500/{crop_id}/*labels.tif",
            f"results_batos/**/{scan}*{crop_id}*/**/*batos*min500*labels.tif",
            f"results_batos/**/{crop_id}*/**/*batos*min500*labels.tif",
        ]
        batos = find_one(batos_patterns)

        # BA-TOS raw / sem size filter
        raw_patterns = [
            f"{base}/{crop_id}_x{x}_y{y}_z{z}_s200/batos_pipeline/05_watershed/*labels.tif",
            f"{base}/{crop_id}_x{x}_y{y}_z{z}_s200/batos_result/05_watershed/*labels.tif",
            f"results_batos/**/{scan}*{crop_id}*/**/batos_pipeline/05_watershed/*labels.tif",
            f"results_batos/**/{crop_id}*/**/batos_pipeline/05_watershed/*labels.tif",
            f"results_batos/**/{crop_id}*/**/*batos_distance_labels.tif",
        ]
        batos_raw = find_one(raw_patterns)

        # Otsu
        otsu_patterns = [
            f"results_batos/baselines/otsu_selected_crops/{crop_id}/*otsu*labels.tif",
            f"results_batos/{scan}_selected_crops/{crop_id}/otsu*/**/*labels.tif",
            f"results_batos/**/{scan}*{crop_id}*/**/*otsu*labels.tif",
            f"results_batos/**/{crop_id}*/**/*otsu*labels.tif",
        ]
        otsu = find_one(otsu_patterns)

        # Min-tree-core+WS
        mintree_patterns = [
            f"results_batos/baselines/macke_mintree_selected_crops/{crop_id}/**/*labels.tif",
            f"results_batos/{scan}_selected_crops/{crop_id}/macke_mintree*/**/*labels.tif",
            f"results_batos/**/{scan}*{crop_id}*/**/*mintree*labels.tif",
            f"results_batos/**/{crop_id}*/**/*mintree*labels.tif",
            f"results_batos/**/{crop_id}*/**/*macke*labels.tif",
        ]
        mintree = find_one(mintree_patterns)

        for method, pred in [
            ("BA-TOS", batos),
            ("BA-TOS raw/no-size-filter", batos_raw),
            ("Otsu+WS", otsu),
            ("Min-tree-core+WS", mintree),
        ]:
            rows.append({
                "scan": scan,
                "crop_id": crop_id,
                "gt": str(gt) if gt else "",
                "method": method,
                "pred": str(pred) if pred else "",
                "gt_found": int(gt is not None),
                "pred_found": int(pred is not None),
            })

out_csv = OUT / "revision_manifest.csv"
with out_csv.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["scan", "crop_id", "gt", "method", "pred", "gt_found", "pred_found"])
    w.writeheader()
    w.writerows(rows)

print(f"wrote: {out_csv}")
print()
print("=== missing entries ===")
missing = [r for r in rows if r["gt_found"] == 0 or r["pred_found"] == 0]
if not missing:
    print("none")
else:
    for r in missing:
        print(r["scan"], r["crop_id"], r["method"], "gt=", r["gt_found"], "pred=", r["pred_found"])
