#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys
import json
import numpy as np
import pandas as pd
import tifffile as tiff

ROOT = Path(".")
SCRIPT_DIR = ROOT / "code/segmentation/batos"

TILE_DIR = ROOT / "results_batos/full_image_pilot_tiles/EFRGP01_01_tile_center_256_halo_64"
PRED_DIR = TILE_DIR / "batos_ls_auto_prediction"

GT_CENTER = PRED_DIR / "07_center_crop/gt_label_center_256.tif"
PRED_CENTER = PRED_DIR / "07_center_crop/batos_ls_auto_pred_center_256.tif"
BATOS_CENTER = PRED_DIR / "07_center_crop/batos_pred_center_256.tif"

METRICS = SCRIPT_DIR / "batos_instance_metrics.py"

OUT_METRICS = PRED_DIR / "08_metrics"
OUT_PREVIEW = PRED_DIR / "09_gmic_preview"
OUT_METRICS.mkdir(parents=True, exist_ok=True)
OUT_PREVIEW.mkdir(parents=True, exist_ok=True)

THRESHOLDS = ["0.25", "0.50", "0.75"]


def require(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def run_metric(gt, pred, threshold):
    tag = threshold.replace(".", "p")
    name = f"pilot_tile_center_256_BA_TOS_LS_AUTO_iou_{tag}"
    out_dir = OUT_METRICS / f"iou_{tag}"
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, str(METRICS),
        "--gt", str(gt),
        "--pred", str(pred),
        "--out-dir", str(out_dir),
        "--min-gt-area", "1000",
        "--ignore-gt-border",
        "--thresholds", threshold,
        "--name", name,
    ]

    print("RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    csv_path = out_dir / f"{name}_instance_metrics.csv"

    if not csv_path.exists():
        raise FileNotFoundError(csv_path)

    row = pd.read_csv(csv_path).iloc[0].to_dict()
    row["threshold"] = float(threshold)
    row["csv"] = str(csv_path)
    return row


def colorize_labels(labels):
    labels = labels.astype(np.int64)
    rgb = np.zeros(labels.shape + (3,), dtype=np.uint8)

    ids = np.unique(labels)
    ids = ids[ids > 0]

    for lab in ids:
        r = int((lab * 37) % 255)
        g = int((lab * 91) % 255)
        b = int((lab * 53) % 255)

        r = max(r, 60)
        g = max(g, 60)
        b = max(b, 60)

        rgb[labels == lab] = (r, g, b)

    return rgb


def center_crop(arr, size=160):
    z, y, x = arr.shape
    z0 = max(0, z // 2 - size // 2)
    y0 = max(0, y // 2 - size // 2)
    x0 = max(0, x // 2 - size // 2)
    return arr[z0:z0+size, y0:y0+size, x0:x0+size]


def make_preview():
    gt = tiff.imread(str(GT_CENTER))
    pred = tiff.imread(str(PRED_CENTER))
    batos = tiff.imread(str(BATOS_CENTER))

    if gt.shape != pred.shape:
        raise RuntimeError(f"Shape mismatch: gt={gt.shape}, pred={pred.shape}")

    gt160 = center_crop(gt, 160)
    batos160 = center_crop(batos, 160)
    pred160 = center_crop(pred, 160)

    gt_rgb = colorize_labels(gt160)
    batos_rgb = colorize_labels(batos160)
    pred_rgb = colorize_labels(pred160)

    tiff.imwrite(str(OUT_PREVIEW / "gt_center_160_rgb_stack.tif"), gt_rgb, photometric="rgb")
    tiff.imwrite(str(OUT_PREVIEW / "batos_center_160_rgb_stack.tif"), batos_rgb, photometric="rgb")
    tiff.imwrite(str(OUT_PREVIEW / "batos_ls_auto_center_160_rgb_stack.tif"), pred_rgb, photometric="rgb")

    gap = np.ones((gt_rgb.shape[0], gt_rgb.shape[1], 8, 3), dtype=np.uint8) * 255
    side = np.concatenate([gt_rgb, gap, batos_rgb, gap, pred_rgb], axis=2)

    side_path = OUT_PREVIEW / "gt_vs_batos_vs_batos_ls_auto_center_160_rgb_stack.tif"
    tiff.imwrite(str(side_path), side, photometric="rgb")

    print()
    print("G'MIC preview written:")
    print(side_path)
    print()
    print("Open with:")
    print(f"gmic {side_path} -append z -display")


def main():
    require(GT_CENTER)
    require(PRED_CENTER)
    require(BATOS_CENTER)
    require(METRICS)

    rows = []

    for th in THRESHOLDS:
        rows.append(run_metric(GT_CENTER, PRED_CENTER, th))

    df = pd.DataFrame(rows)
    crop_metrics_csv = OUT_METRICS / "pilot_tile_center_256_BA_TOS_LS_AUTO_multi_iou_metrics.csv"
    df.to_csv(crop_metrics_csv, index=False)

    row025 = df[df["threshold"] == 0.25].iloc[0]
    row050 = df[df["threshold"] == 0.50].iloc[0]
    row075 = df[df["threshold"] == 0.75].iloc[0]

    summary = pd.DataFrame([{
        "method": "BA-TOS-LS-AUTO",
        "tile": "EFRGP01_01_center_256_halo_64",
        "P@0.50": row050["precision"],
        "R@0.50": row050["recall"],
        "F1@0.25": row025["f1"],
        "F1@0.50": row050["f1"],
        "F1@0.75": row075["f1"],
        "overseg": row050["overseg_ratio_pred_eval_gt"],
        "tp@0.50": int(row050["tp"]),
        "fp@0.50": int(row050["fp"]),
        "fn@0.50": int(row050["fn"]),
    }])

    summary_csv = OUT_METRICS / "pilot_tile_center_256_BA_TOS_LS_AUTO_summary.csv"
    summary.to_csv(summary_csv, index=False)

    print()
    print("=== Pilot tile BA-TOS-LS-AUTO summary ===")
    print(summary.to_string(index=False))
    print()
    print("wrote:", crop_metrics_csv)
    print("wrote:", summary_csv)

    make_preview()

    manifest = {
        "gt_center": str(GT_CENTER),
        "batos_center": str(BATOS_CENTER),
        "batos_ls_auto_center": str(PRED_CENTER),
        "multi_iou_metrics": str(crop_metrics_csv),
        "summary": str(summary_csv),
        "gmic_preview": str(OUT_PREVIEW / "gt_vs_batos_vs_batos_ls_auto_center_160_rgb_stack.tif"),
    }

    with open(PRED_DIR / "pilot_tile_finalize_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    print("wrote:", PRED_DIR / "pilot_tile_finalize_manifest.json")


if __name__ == "__main__":
    main()
