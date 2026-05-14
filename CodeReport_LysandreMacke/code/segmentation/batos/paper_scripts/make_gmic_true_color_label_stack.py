#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import tifffile as tiff

IN_DIR = Path("results_batos/full_image_pilot_tiles/EFRGP01_01_tile_center_256_halo_64")
OUT = IN_DIR / "gmic_color_3d"
OUT.mkdir(parents=True, exist_ok=True)

LABEL = IN_DIR / "label_center_256.tif"
PREVIEW_SIZE = 160


def center_crop(arr, size):
    z, y, x = arr.shape
    z0 = max(0, z // 2 - size // 2)
    y0 = max(0, y // 2 - size // 2)
    x0 = max(0, x // 2 - size // 2)
    return arr[z0:z0+size, y0:y0+size, x0:x0+size]


def label_to_rgb_stack(labels):
    labels = labels.astype(np.int64)
    z, y, x = labels.shape
    rgb = np.zeros((z, y, x, 3), dtype=np.uint8)

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


def compute_surface(labels):
    labels = labels.astype(np.int64)
    surf = np.zeros(labels.shape, dtype=bool)

    surf[1:, :, :] |= labels[1:, :, :] != labels[:-1, :, :]
    surf[:-1, :, :] |= labels[1:, :, :] != labels[:-1, :, :]

    surf[:, 1:, :] |= labels[:, 1:, :] != labels[:, :-1, :]
    surf[:, :-1, :] |= labels[:, 1:, :] != labels[:, :-1, :]

    surf[:, :, 1:] |= labels[:, :, 1:] != labels[:, :, :-1]
    surf[:, :, :-1] |= labels[:, :, 1:] != labels[:, :, :-1]

    surf &= labels > 0
    return surf


def make_surface_rgb_stack(labels):
    full_rgb = label_to_rgb_stack(labels)
    surf = compute_surface(labels)
    out = np.zeros_like(full_rgb)
    out[surf] = full_rgb[surf]
    return out


def main():
    labels = tiff.imread(str(LABEL))
    labels = center_crop(labels, PREVIEW_SIZE)

    rgb_full = label_to_rgb_stack(labels)
    rgb_surface = make_surface_rgb_stack(labels)

    tiff.imwrite(
        str(OUT / "label_center_160_rgb_stack.tif"),
        rgb_full,
        photometric="rgb"
    )

    tiff.imwrite(
        str(OUT / "label_center_160_surface_rgb_stack.tif"),
        rgb_surface,
        photometric="rgb"
    )

    print("wrote:", OUT / "label_center_160_rgb_stack.tif")
    print("wrote:", OUT / "label_center_160_surface_rgb_stack.tif")
    print("shape full   :", rgb_full.shape)
    print("shape surface:", rgb_surface.shape)


if __name__ == "__main__":
    main()
