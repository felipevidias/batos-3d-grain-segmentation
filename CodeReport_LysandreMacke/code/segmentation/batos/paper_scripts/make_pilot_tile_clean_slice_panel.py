#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt

ROOT = Path(".")
TILE = ROOT / "results_batos/full_image_pilot_tiles/EFRGP01_01_tile_center_256_halo_64"

GRAY = TILE / "gray_center_256.tif"
GT = TILE / "label_center_256.tif"

BATOS = TILE / "batos_ls_auto_tile_rerun_maxmarkers_1800/07_center_crop/batos_pred_center_256.tif"
AUTO = TILE / "batos_ls_auto_tile_rerun_maxmarkers_1800/07_center_crop/batos_ls_auto_tile_pred_center_256.tif"

OUT = TILE / "batos_ls_auto_tile_rerun_maxmarkers_1800/10_slice_panels"
OUT.mkdir(parents=True, exist_ok=True)


def read(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D, got {arr.shape}: {path}")
    return arr


def normalize_gray(img):
    img = img.astype(np.float32)
    p1, p99 = np.percentile(img, [1, 99])
    img = (img - p1) / (p99 - p1 + 1e-8)
    return np.clip(img, 0, 1)


def colorize_labels_2d(labels):
    labels = labels.astype(np.int64)
    rgb = np.zeros(labels.shape + (3,), dtype=np.float32)

    ids = np.unique(labels)
    ids = ids[ids > 0]

    for lab in ids:
        r = ((lab * 37) % 255) / 255.0
        g = ((lab * 91) % 255) / 255.0
        b = ((lab * 53) % 255) / 255.0

        r = max(r, 0.25)
        g = max(g, 0.25)
        b = max(b, 0.25)

        rgb[labels == lab] = (r, g, b)

    return rgb


def overlay(gray, labels, alpha=0.55):
    g = normalize_gray(gray)
    base = np.repeat(g[..., None], 3, axis=-1)
    col = colorize_labels_2d(labels)
    mask = labels > 0
    out = base.copy()
    out[mask] = (1 - alpha) * base[mask] + alpha * col[mask]
    return np.clip(out, 0, 1)


def choose_slices(gt, pred, n=6):
    # escolhe slices com bastante GT e diferença entre GT/pred
    scores = []
    for z in range(gt.shape[0]):
        gt_count = np.count_nonzero(gt[z])
        pred_count = np.count_nonzero(pred[z])
        diff = abs(gt_count - pred_count)
        scores.append((gt_count + 0.5 * diff, z))
    scores = sorted(scores, reverse=True)
    zs = [z for _, z in scores[:n]]
    return sorted(zs)


def make_panel(z, gray, gt, batos, auto):
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))

    axes[0].imshow(normalize_gray(gray[z]), cmap="gray")
    axes[0].set_title(f"Grayscale z={z}")

    axes[1].imshow(colorize_labels_2d(gt[z]))
    axes[1].set_title("Reference")

    axes[2].imshow(overlay(gray[z], batos[z]))
    axes[2].set_title("BA-TOS")

    axes[3].imshow(overlay(gray[z], auto[z]))
    axes[3].set_title("BA-TOS-LS-AUTO-TILE")

    for ax in axes:
        ax.axis("off")

    fig.tight_layout()
    out = OUT / f"pilot_tile_panel_z{z:03d}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote:", out)


def main():
    gray = read(GRAY)
    gt = read(GT)
    batos = read(BATOS)
    auto = read(AUTO)

    print("gray:", gray.shape)
    print("gt:", gt.shape, "labels:", len(np.unique(gt)) - 1)
    print("batos:", batos.shape, "labels:", len(np.unique(batos)) - 1)
    print("auto:", auto.shape, "labels:", len(np.unique(auto)) - 1)

    zs = choose_slices(gt, auto, n=8)
    print("selected slices:", zs)

    for z in zs:
        make_panel(z, gray, gt, batos, auto)


if __name__ == "__main__":
    main()
