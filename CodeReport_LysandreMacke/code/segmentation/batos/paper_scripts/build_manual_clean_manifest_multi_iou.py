#!/usr/bin/env python3
from pathlib import Path
import csv

OUT = Path("results_batos/paper_tables/final_revision_clean")
OUT.mkdir(parents=True, exist_ok=True)

CROPS = {
    "crop02": ("600", "800", "300"),
    "crop03": ("900", "600", "300"),
    "crop04": ("500", "500", "700"),
    "crop05": ("800", "900", "700"),
    "crop06": ("350", "650", "500"),
    "crop07": ("700", "400", "500"),
}

rows = []

def add(scan, crop_id, gt, method, pred):
    gt = Path(gt)
    pred = Path(pred)
    rows.append({
        "scan": scan,
        "crop_id": crop_id,
        "method": method,
        "gt": str(gt),
        "pred": str(pred),
        "gt_found": int(gt.exists()),
        "pred_found": int(pred.exists()),
    })

for crop_id, (x, y, z) in CROPS.items():

    # =========================
    # EFRGP01_00
    # =========================
    if crop_id in ["crop02", "crop03", "crop04", "crop05"]:
        base00 = Path(f"results_batos/EFRGP01_00_multicrop/{crop_id}_x{x}_y{y}_z{z}_s200")
    else:
        base00 = Path(f"results_batos/EFRGP01_00_extra_crops/{crop_id}_x{x}_y{y}_z{z}_s200")

    gt00 = base00 / "00_input" / f"{crop_id}_x{x}_y{y}_z{z}_s200_label_crop.tif"

    batos00 = Path(f"results_batos/paper_tables/batos_size_filter_sweep/min500/{crop_id}/{crop_id}_batos_min500_labels.tif")
    otsu00 = Path(f"results_batos/baselines/otsu_selected_crops/{crop_id}/{crop_id}_otsu_ws_labels.tif")
    macke00 = Path(f"results_batos/baselines/macke_mintree_selected_crops/{crop_id}/02_watershed/macke_mintree_result.tif")

    add("EFRGP01_00", crop_id, gt00, "BA-TOS", batos00)
    add("EFRGP01_00", crop_id, gt00, "Otsu+WS", otsu00)
    add("EFRGP01_00", crop_id, gt00, "Min-tree-core+WS", macke00)

    # =========================
    # EFRGP01_01
    # =========================
    base01 = Path(f"results_batos/EFRGP01_01_selected_crops/{crop_id}_x{x}_y{y}_z{z}_s200")

    gt01 = base01 / "00_input" / f"EFRGP01_01_{crop_id}_x{x}_y{y}_z{z}_s200_label_crop.tif"

    batos01 = base01 / "batos_pipeline" / "05_watershed" / "batos_distance_filtered_min500_labels.tif"
    otsu01 = base01 / "otsu_ws" / f"{crop_id}_otsu_ws_labels.tif"
    macke01 = base01 / "macke_mintree" / "02_watershed" / "macke_mintree_result.tif"

    add("EFRGP01_01", crop_id, gt01, "BA-TOS", batos01)
    add("EFRGP01_01", crop_id, gt01, "Otsu+WS", otsu01)
    add("EFRGP01_01", crop_id, gt01, "Min-tree-core+WS", macke01)

out_csv = OUT / "manual_clean_manifest_multi_iou.csv"
with out_csv.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=["scan", "crop_id", "method", "gt", "pred", "gt_found", "pred_found"]
    )
    w.writeheader()
    w.writerows(rows)

print("wrote:", out_csv)

print("\nCounts by method:")
for method in ["Min-tree-core+WS", "Otsu+WS", "BA-TOS"]:
    sub = [r for r in rows if r["method"] == method]
    ok = [r for r in sub if r["gt_found"] == 1 and r["pred_found"] == 1]
    print(f"{method}: {len(ok)}/{len(sub)}")

print("\nMissing:")
missing = [r for r in rows if r["gt_found"] == 0 or r["pred_found"] == 0]
if not missing:
    print("none")
else:
    for r in missing:
        print(r["scan"], r["crop_id"], r["method"])
        print("  gt_found  :", r["gt_found"], r["gt"])
        print("  pred_found:", r["pred_found"], r["pred"])
