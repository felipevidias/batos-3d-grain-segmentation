#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import tifffile as tiff
import pandas as pd

ROOT = Path(".")
OUT = ROOT / "results_batos/paper_figures/batos_ls_color_qc"
OUT.mkdir(parents=True, exist_ok=True)

CASES = [
    {
        "name": "EFRGP01_01_crop04",
        "gray": "results_batos/EFRGP01_01_selected_crops/crop04_x500_y500_z700_s200/00_input/EFRGP01_01_crop04_x500_y500_z700_s200_gray_crop_8u.tif",
        "gt": "results_batos/EFRGP01_01_selected_crops/crop04_x500_y500_z700_s200/00_input/EFRGP01_01_crop04_x500_y500_z700_s200_label_crop.tif",
        "orig": "results_batos/EFRGP01_01_selected_crops/crop04_x500_y500_z700_s200/batos_pipeline/05_watershed/batos_distance_filtered_min500_labels.tif",
        "ls": "results_batos/paper_tables/batos_local_split_large_labels/EFRGP01_01/crop04/EFRGP01_01_crop04_batos_local_split_labels.tif",
    },
    {
        "name": "EFRGP01_01_crop07",
        "gray": "results_batos/EFRGP01_01_selected_crops/crop07_x700_y400_z500_s200/00_input/EFRGP01_01_crop07_x700_y400_z500_s200_gray_crop_8u.tif",
        "gt": "results_batos/EFRGP01_01_selected_crops/crop07_x700_y400_z500_s200/00_input/EFRGP01_01_crop07_x700_y400_z500_s200_label_crop.tif",
        "orig": "results_batos/EFRGP01_01_selected_crops/crop07_x700_y400_z500_s200/batos_pipeline/05_watershed/batos_distance_filtered_min500_labels.tif",
        "ls": "results_batos/paper_tables/batos_local_split_large_labels/EFRGP01_01/crop07/EFRGP01_01_crop07_batos_local_split_labels.tif",
    },
    {
        "name": "EFRGP01_00_crop05",
        "gray": "results_batos/EFRGP01_00_multicrop/crop05_x800_y900_z700_s200/00_input/crop05_x800_y900_z700_s200_gray_crop_8u.tif",
        "gt": "results_batos/EFRGP01_00_multicrop/crop05_x800_y900_z700_s200/00_input/crop05_x800_y900_z700_s200_label_crop.tif",
        "orig": "results_batos/paper_tables/batos_size_filter_sweep/min500/crop05/crop05_batos_min500_labels.tif",
        "ls": "results_batos/paper_tables/batos_local_split_large_labels/EFRGP01_00/crop05/EFRGP01_00_crop05_batos_local_split_labels.tif",
    },
]


def read(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    arr = tiff.imread(str(path))
    return arr


def make_palette(max_label, seed=12345):
    rng = np.random.default_rng(seed)
    palette = np.zeros((max_label + 1, 3), dtype=np.uint8)
    if max_label > 0:
        palette[1:] = rng.integers(40, 255, size=(max_label, 3), dtype=np.uint8)
    palette[0] = [0, 0, 0]
    return palette


def labels_to_rgb(labels, seed=12345):
    labels = labels.astype(np.int64, copy=False)
    max_label = int(labels.max())
    palette = make_palette(max_label, seed=seed)
    rgb = palette[labels]
    return rgb.astype(np.uint8)


def normalize_gray(gray):
    gray = gray.astype(np.float32)
    lo, hi = np.percentile(gray, [1, 99])
    if hi <= lo:
        hi = gray.max() if gray.max() > lo else lo + 1
    gray = np.clip((gray - lo) / (hi - lo), 0, 1)
    return (gray * 255).astype(np.uint8)


def overlay_labels(gray, labels, seed=12345, alpha=0.45):
    g = normalize_gray(gray)
    base = np.stack([g, g, g], axis=-1).astype(np.float32)
    color = labels_to_rgb(labels, seed=seed).astype(np.float32)

    mask = labels > 0
    out = base.copy()
    out[mask] = (1 - alpha) * base[mask] + alpha * color[mask]
    return np.clip(out, 0, 255).astype(np.uint8)


def save_mip_png(volume_rgb, out_png, axis=0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if volume_rgb.ndim == 4:
        mip = volume_rgb.max(axis=axis)
    else:
        mip = volume_rgb.max(axis=axis)
        mip = np.stack([mip, mip, mip], axis=-1)

    plt.figure(figsize=(6, 6))
    plt.imshow(mip)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def main():
    rows = []

    for case in CASES:
        name = case["name"]
        case_out = OUT / name
        case_out.mkdir(parents=True, exist_ok=True)

        print("[case]", name)

        gray = read(case["gray"])
        gt = read(case["gt"])
        orig = read(case["orig"])
        ls = read(case["ls"])

        gt_rgb = labels_to_rgb(gt, seed=10)
        orig_rgb = labels_to_rgb(orig, seed=20)
        ls_rgb = labels_to_rgb(ls, seed=30)

        gt_ov = overlay_labels(gray, gt, seed=10, alpha=0.45)
        orig_ov = overlay_labels(gray, orig, seed=20, alpha=0.45)
        ls_ov = overlay_labels(gray, ls, seed=30, alpha=0.45)

        files = {
            "gray": case_out / f"{name}_gray.tif",
            "gt_rgb": case_out / f"{name}_gt_color_rgb.tif",
            "batos_rgb": case_out / f"{name}_batos_color_rgb.tif",
            "batos_ls_rgb": case_out / f"{name}_batos_ls_color_rgb.tif",
            "gt_overlay": case_out / f"{name}_gt_overlay_rgb.tif",
            "batos_overlay": case_out / f"{name}_batos_overlay_rgb.tif",
            "batos_ls_overlay": case_out / f"{name}_batos_ls_overlay_rgb.tif",
        }

        tiff.imwrite(files["gray"], gray.astype(np.uint8))
        tiff.imwrite(files["gt_rgb"], gt_rgb)
        tiff.imwrite(files["batos_rgb"], orig_rgb)
        tiff.imwrite(files["batos_ls_rgb"], ls_rgb)
        tiff.imwrite(files["gt_overlay"], gt_ov)
        tiff.imwrite(files["batos_overlay"], orig_ov)
        tiff.imwrite(files["batos_ls_overlay"], ls_ov)

        save_mip_png(gt_ov, case_out / f"{name}_gt_overlay_MIP.png")
        save_mip_png(orig_ov, case_out / f"{name}_batos_overlay_MIP.png")
        save_mip_png(ls_ov, case_out / f"{name}_batos_ls_overlay_MIP.png")

        rows.append({
            "case": name,
            **{k: str(v) for k, v in files.items()}
        })

        print("wrote:", case_out)

    pd.DataFrame(rows).to_csv(OUT / "color_qc_manifest.csv", index=False)
    print("wrote:", OUT / "color_qc_manifest.csv")


if __name__ == "__main__":
    main()
