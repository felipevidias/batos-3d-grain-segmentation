#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import tifffile as tiff
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skimage.color import label2rgb

ROOT = Path(".")
OUT = ROOT / "results_batos/paper_figures/batos_ls_png_qc"
OUT.mkdir(parents=True, exist_ok=True)

CASES = [
    {
        "name": "EFRGP01_01_crop04",
        "gray": "results_batos/EFRGP01_01_selected_crops/crop04_x500_y500_z700_s200/00_input/EFRGP01_01_crop04_x500_y500_z700_s200_gray_crop_8u.tif",
        "gt": "results_batos/EFRGP01_01_selected_crops/crop04_x500_y500_z700_s200/00_input/EFRGP01_01_crop04_x500_y500_z700_s200_label_crop.tif",
        "batos": "results_batos/EFRGP01_01_selected_crops/crop04_x500_y500_z700_s200/batos_pipeline/05_watershed/batos_distance_filtered_min500_labels.tif",
        "ls": "results_batos/paper_tables/batos_local_split_large_labels/EFRGP01_01/crop04/EFRGP01_01_crop04_batos_local_split_labels.tif",
    },
    {
        "name": "EFRGP01_01_crop07",
        "gray": "results_batos/EFRGP01_01_selected_crops/crop07_x700_y400_z500_s200/00_input/EFRGP01_01_crop07_x700_y400_z500_s200_gray_crop_8u.tif",
        "gt": "results_batos/EFRGP01_01_selected_crops/crop07_x700_y400_z500_s200/00_input/EFRGP01_01_crop07_x700_y400_z500_s200_label_crop.tif",
        "batos": "results_batos/EFRGP01_01_selected_crops/crop07_x700_y400_z500_s200/batos_pipeline/05_watershed/batos_distance_filtered_min500_labels.tif",
        "ls": "results_batos/paper_tables/batos_local_split_large_labels/EFRGP01_01/crop07/EFRGP01_01_crop07_batos_local_split_labels.tif",
    },
    {
        "name": "EFRGP01_00_crop05",
        "gray": "results_batos/EFRGP01_00_multicrop/crop05_x800_y900_z700_s200/00_input/crop05_x800_y900_z700_s200_gray_crop_8u.tif",
        "gt": "results_batos/EFRGP01_00_multicrop/crop05_x800_y900_z700_s200/00_input/crop05_x800_y900_z700_s200_label_crop.tif",
        "batos": "results_batos/paper_tables/batos_size_filter_sweep/min500/crop05/crop05_batos_min500_labels.tif",
        "ls": "results_batos/paper_tables/batos_local_split_large_labels/EFRGP01_00/crop05/EFRGP01_00_crop05_batos_local_split_labels.tif",
    },
]


def read(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return tiff.imread(str(path))


def norm01(img):
    img = img.astype(np.float32)
    lo, hi = np.percentile(img, [1, 99])
    if hi <= lo:
        hi = img.max() if img.max() > lo else lo + 1
    return np.clip((img - lo) / (hi - lo), 0, 1)


def labels_overlay(labels2d, gray2d):
    g = norm01(gray2d)
    rgb = label2rgb(
        labels2d.astype(np.int32),
        image=g,
        bg_label=0,
        alpha=0.55,
        saturation=0,
        kind="overlay",
    )
    return np.clip(rgb, 0, 1)


def choose_slices(gt, batos, ls, n=5):
    # Escolhe fatias com maior quantidade de diferença entre BA-TOS e BA-TOS-LS,
    # e também com bastante GT para ser visível.
    diff = (batos != ls)
    scores = []
    for z in range(gt.shape[0]):
        score = int(diff[z].sum()) + 0.05 * int((gt[z] > 0).sum())
        scores.append(score)
    scores = np.array(scores)
    idx = np.argsort(scores)[::-1]
    chosen = []
    for z in idx:
        z = int(z)
        if len(chosen) == 0 or all(abs(z - c) >= 8 for c in chosen):
            chosen.append(z)
        if len(chosen) >= n:
            break
    return sorted(chosen)


def make_panel(case):
    name = case["name"]
    outdir = OUT / name
    outdir.mkdir(parents=True, exist_ok=True)

    gray = read(case["gray"])
    gt = read(case["gt"])
    batos = read(case["batos"])
    ls = read(case["ls"])

    slices = choose_slices(gt, batos, ls, n=6)
    print(name, "selected slices:", slices)

    for z in slices:
        gray2 = gray[z]
        gt2 = gt[z]
        batos2 = batos[z]
        ls2 = ls[z]

        # Diferença: vermelho = BA-TOS original, verde = BA-TOS-LS, amarelo = comum
        bmask = batos2 > 0
        lmask = ls2 > 0
        g = norm01(gray2)
        diff_rgb = np.stack([g, g, g], axis=-1)
        only_b = bmask & (~lmask)
        only_l = lmask & (~bmask)
        both = bmask & lmask
        diff_rgb[only_b] = [1, 0, 0]
        diff_rgb[only_l] = [0, 1, 0]
        diff_rgb[both] = [1, 1, 0]

        fig, axes = plt.subplots(1, 5, figsize=(18, 4))
        axes[0].imshow(g, cmap="gray")
        axes[0].set_title("gray")

        axes[1].imshow(labels_overlay(gt2, gray2))
        axes[1].set_title("GT overlay")

        axes[2].imshow(labels_overlay(batos2, gray2))
        axes[2].set_title("BA-TOS overlay")

        axes[3].imshow(labels_overlay(ls2, gray2))
        axes[3].set_title("BA-TOS-LS overlay")

        axes[4].imshow(diff_rgb)
        axes[4].set_title("diff: red=B, green=LS, yellow=both")

        for ax in axes:
            ax.axis("off")

        fig.suptitle(f"{name} | z={z}", fontsize=11)
        fig.tight_layout()

        outpng = outdir / f"{name}_z{z:03d}_qc_panel.png"
        fig.savefig(outpng, dpi=180)
        plt.close(fig)
        print("wrote:", outpng)

    # Também salva um MIP simples só para prévia
    gt_mip = (gt > 0).max(axis=0).astype(np.uint8)
    b_mip = (batos > 0).max(axis=0).astype(np.uint8)
    l_mip = (ls > 0).max(axis=0).astype(np.uint8)
    g_mip = gray.max(axis=0)

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    axes[0].imshow(norm01(g_mip), cmap="gray")
    axes[0].set_title("gray MIP")
    axes[1].imshow(gt_mip, cmap="gray")
    axes[1].set_title("GT MIP mask")
    axes[2].imshow(b_mip, cmap="gray")
    axes[2].set_title("BA-TOS MIP mask")
    axes[3].imshow(l_mip, cmap="gray")
    axes[3].set_title("BA-TOS-LS MIP mask")
    for ax in axes:
        ax.axis("off")
    fig.tight_layout()
    outpng = outdir / f"{name}_mip_masks.png"
    fig.savefig(outpng, dpi=180)
    plt.close(fig)


def main():
    for case in CASES:
        make_panel(case)
    print("done:", OUT)


if __name__ == "__main__":
    main()
