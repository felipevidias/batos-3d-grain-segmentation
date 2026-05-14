#!/usr/bin/env python3
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import tifffile as tiff

ROOT = Path(".")
MANIFEST = ROOT / "results_batos/paper_tables/batos_local_split_auto/local_split_auto_manifest.csv"

OUT_ROOT = ROOT / "results_batos/fiji_3d_view_batos_ls_auto"
OUT_ROOT.mkdir(parents=True, exist_ok=True)


def read_label(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise ValueError(f"Expected 3D scalar label, got {arr.shape}: {path}")
    return arr.astype(np.int32)


def label_to_rgb_volume(labels):
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


def process_one(row):
    scan = str(row["scan"])
    crop_id = str(row["crop_id"])
    name = f"{scan}_{crop_id}"

    gt_path = Path(str(row["gt"]))
    batos_path = Path(str(row["original_pred"]))
    auto_path = Path(str(row["split_pred"]))

    out = OUT_ROOT / name
    out.mkdir(parents=True, exist_ok=True)

    # Escalares reais
    shutil.copy(gt_path, out / "01_GT_labels_3D_scalar.tif")
    shutil.copy(batos_path, out / "02_BA_TOS_labels_3D_scalar.tif")
    shutil.copy(auto_path, out / "03_BA_TOS_LS_AUTO_labels_3D_scalar.tif")

    gt = read_label(gt_path)
    batos = read_label(batos_path)
    auto = read_label(auto_path)

    gt_rgb = label_to_rgb_volume(gt)
    batos_rgb = label_to_rgb_volume(batos)
    auto_rgb = label_to_rgb_volume(auto)

    tiff.imwrite(str(out / "01_GT_labels_3D_RGB.tif"), gt_rgb, photometric="rgb")
    tiff.imwrite(str(out / "02_BA_TOS_labels_3D_RGB.tif"), batos_rgb, photometric="rgb")
    tiff.imwrite(str(out / "03_BA_TOS_LS_AUTO_labels_3D_RGB.tif"), auto_rgb, photometric="rgb")

    z, y, x, c = gt_rgb.shape
    gap = np.ones((z, y, 8, c), dtype=np.uint8) * 255

    side = np.concatenate([gt_rgb, gap, batos_rgb, gap, auto_rgb], axis=2)

    side_out = out / "00_GT_BA_TOS_BA_TOS_LS_AUTO_side_by_side_3D_RGB.tif"
    tiff.imwrite(str(side_out), side, photometric="rgb")

    print("[OK]", side_out)


def main():
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Missing manifest: {MANIFEST}")

    df = pd.read_csv(MANIFEST)

    for _, row in df.iterrows():
        process_one(row)

    print()
    print("Fiji volumes written to:", OUT_ROOT)


if __name__ == "__main__":
    main()
