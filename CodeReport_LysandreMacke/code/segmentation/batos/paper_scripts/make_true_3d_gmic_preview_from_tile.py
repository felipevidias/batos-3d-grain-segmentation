#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import tifffile as tiff

IN_DIR = Path("results_batos/full_image_pilot_tiles/EFRGP01_01_tile_center_256_halo_64")
OUT = IN_DIR / "gmic_true_3d_preview"
OUT.mkdir(parents=True, exist_ok=True)

LABEL = IN_DIR / "label_center_256.tif"
GRAY = IN_DIR / "gray_center_256.tif"

PREVIEW_SIZE = 160


def center_crop(arr, size):
    z, y, x = arr.shape
    z0 = max(0, z // 2 - size // 2)
    y0 = max(0, y // 2 - size // 2)
    x0 = max(0, x // 2 - size // 2)
    return arr[z0:z0+size, y0:y0+size, x0:x0+size]


def normalize_gray(arr):
    arr = arr.astype(np.float32)
    p1, p99 = np.percentile(arr, [1, 99])
    arr = (arr - p1) / (p99 - p1 + 1e-8)
    arr = np.clip(arr, 0, 1)
    return (arr * 255).astype(np.uint8)


def compute_surface(labels):
    labels = labels.astype(np.int64)
    surf = np.zeros(labels.shape, dtype=np.uint8)

    diff = np.zeros(labels.shape, dtype=bool)

    diff[1:, :, :] |= labels[1:, :, :] != labels[:-1, :, :]
    diff[:-1, :, :] |= labels[1:, :, :] != labels[:-1, :, :]

    diff[:, 1:, :] |= labels[:, 1:, :] != labels[:, :-1, :]
    diff[:, :-1, :] |= labels[:, 1:, :] != labels[:, :-1, :]

    diff[:, :, 1:] |= labels[:, :, 1:] != labels[:, :, :-1]
    diff[:, :, :-1] |= labels[:, :, 1:] != labels[:, :, :-1]

    diff &= labels > 0
    surf[diff] = 255

    return surf


def make_cutaway(volume, z_frac=0.35, y_frac=0.20):
    vol = volume.copy()
    z, y, x = vol.shape

    vol[:int(z*z_frac), :, :] = 0
    vol[:, :int(y*y_frac), :] = 0

    return vol


def main():
    if not LABEL.exists():
        raise FileNotFoundError(LABEL)
    if not GRAY.exists():
        raise FileNotFoundError(GRAY)

    labels = tiff.imread(str(LABEL))
    gray = tiff.imread(str(GRAY))

    if labels.shape != gray.shape:
        raise RuntimeError(f"Shape mismatch: labels={labels.shape}, gray={gray.shape}")

    labels = center_crop(labels, PREVIEW_SIZE)
    gray = center_crop(gray, PREVIEW_SIZE)

    gray_u8 = normalize_gray(gray)

    # volume escalar com ids de labels
    labels_scalar = labels.astype(np.uint16)

    # binário 3D só de superfície
    surface = compute_surface(labels)

    # superfície com corte
    surface_cut = make_cutaway(surface)

    # versão densidade simples do label (binário preenchido)
    filled = (labels > 0).astype(np.uint8) * 255

    tiff.imwrite(str(OUT / "01_gray_center_160_u8_3d.tif"), gray_u8)
    tiff.imwrite(str(OUT / "02_label_center_160_scalar_3d.tif"), labels_scalar)
    tiff.imwrite(str(OUT / "03_label_center_160_filled_binary_3d.tif"), filled)
    tiff.imwrite(str(OUT / "04_label_center_160_surface_binary_3d.tif"), surface)
    tiff.imwrite(str(OUT / "05_label_center_160_surface_cutaway_binary_3d.tif"), surface_cut)

    print("wrote:", OUT / "01_gray_center_160_u8_3d.tif")
    print("wrote:", OUT / "02_label_center_160_scalar_3d.tif")
    print("wrote:", OUT / "03_label_center_160_filled_binary_3d.tif")
    print("wrote:", OUT / "04_label_center_160_surface_binary_3d.tif")
    print("wrote:", OUT / "05_label_center_160_surface_cutaway_binary_3d.tif")

    print()
    print("shapes:")
    print("gray   :", gray_u8.shape, gray_u8.dtype)
    print("labels :", labels_scalar.shape, labels_scalar.dtype)
    print("filled :", filled.shape, filled.dtype)
    print("surface:", surface.shape, surface.dtype)
    print("cut    :", surface_cut.shape, surface_cut.dtype)


if __name__ == "__main__":
    main()
