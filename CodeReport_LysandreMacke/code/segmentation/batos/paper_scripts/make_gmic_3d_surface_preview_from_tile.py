#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import tifffile as tiff


IN_DIR = Path("results_batos/full_image_pilot_tiles/EFRGP01_01_tile_center_256_halo_64")
OUT = IN_DIR / "gmic_3d_preview"
OUT.mkdir(parents=True, exist_ok=True)

GRAY = IN_DIR / "gray_center_256.tif"
LABEL = IN_DIR / "label_center_256.tif"

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


def label_to_rgb(labels, mask=None):
    labels = labels.astype(np.int64)
    rgb = np.zeros(labels.shape + (3,), dtype=np.uint8)

    ids = np.unique(labels)
    ids = ids[ids > 0]

    for lab in ids:
        r = int((lab * 37) % 255)
        g = int((lab * 91) % 255)
        b = int((lab * 53) % 255)

        r = max(r, 70)
        g = max(g, 70)
        b = max(b, 70)

        lab_mask = labels == lab

        if mask is not None:
            lab_mask = lab_mask & mask

        rgb[lab_mask] = (r, g, b)

    return rgb


def compute_label_surface(labels):
    """
    Marca voxels de superfície/fronteira dos labels.
    Isso evita visualizar o volume preenchido inteiro.
    """
    labels = labels.astype(np.int64)
    surf = np.zeros(labels.shape, dtype=bool)

    # diferenças em z
    surf[1:, :, :] |= labels[1:, :, :] != labels[:-1, :, :]
    surf[:-1, :, :] |= labels[1:, :, :] != labels[:-1, :, :]

    # diferenças em y
    surf[:, 1:, :] |= labels[:, 1:, :] != labels[:, :-1, :]
    surf[:, :-1, :] |= labels[:, 1:, :] != labels[:, :-1, :]

    # diferenças em x
    surf[:, :, 1:] |= labels[:, :, 1:] != labels[:, :, :-1]
    surf[:, :, :-1] |= labels[:, :, 1:] != labels[:, :, :-1]

    surf &= labels > 0
    return surf


def make_cutaway(rgb):
    """
    Remove uma parte frontal para enxergar o interior do bloco.
    """
    cut = rgb.copy()
    z, y, x, c = cut.shape

    # abre uma janela removendo uma fatia frontal/lateral
    cut[: int(0.35 * z), :, :, :] = 0
    cut[:, : int(0.20 * y), :, :] = 0

    return cut


def main():
    if not GRAY.exists():
        raise FileNotFoundError(GRAY)

    if not LABEL.exists():
        raise FileNotFoundError(LABEL)

    gray = tiff.imread(str(GRAY))
    labels = tiff.imread(str(LABEL))

    if gray.shape != labels.shape:
        raise RuntimeError(f"Shape mismatch: gray={gray.shape}, labels={labels.shape}")

    gray = center_crop(gray, PREVIEW_SIZE)
    labels = center_crop(labels, PREVIEW_SIZE)

    print("preview shape:", labels.shape)
    print("labels:", int(len(np.unique(labels)) - 1))
    print("nonzero voxels:", int(np.count_nonzero(labels)))

    gray_u8 = normalize_gray(gray)

    surface = compute_label_surface(labels)

    solid_rgb = label_to_rgb(labels)
    surface_rgb = label_to_rgb(labels, mask=surface)
    surface_cutaway_rgb = make_cutaway(surface_rgb)

    tiff.imwrite(str(OUT / "01_gray_center_160_u8.tif"), gray_u8)
    tiff.imwrite(str(OUT / "02_label_center_160_scalar.tif"), labels.astype(np.uint16))
    tiff.imwrite(str(OUT / "03_label_center_160_solid_rgb.tif"), solid_rgb, photometric="rgb")
    tiff.imwrite(str(OUT / "04_label_center_160_surface_rgb.tif"), surface_rgb, photometric="rgb")
    tiff.imwrite(str(OUT / "05_label_center_160_surface_cutaway_rgb.tif"), surface_cutaway_rgb, photometric="rgb")

    print()
    print("Arquivos para G'MIC:")
    print(OUT / "01_gray_center_160_u8.tif")
    print(OUT / "02_label_center_160_scalar.tif")
    print(OUT / "03_label_center_160_solid_rgb.tif")
    print(OUT / "04_label_center_160_surface_rgb.tif")
    print(OUT / "05_label_center_160_surface_cutaway_rgb.tif")


if __name__ == "__main__":
    main()
