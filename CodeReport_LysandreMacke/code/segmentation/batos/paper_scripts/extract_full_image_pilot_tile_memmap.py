#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import tifffile as tiff

GRAY = Path("/mnt/win_p3/Users/USER/Downloads/PhDUpload/PhDUpload/EFR/Greyscale/EFRGP01_01.tif")
LABEL = Path("/mnt/win_p3/Users/USER/Downloads/PhDUpload/PhDUpload/EFR/Label/extracted_EFRGP01_01/EFRGP01_01.tif")

OUT = Path("results_batos/full_image_pilot_tiles/EFRGP01_01_tile_center_256_halo_64")
OUT.mkdir(parents=True, exist_ok=True)

CENTER_SIZE = 256
HALO = 64
PROC_SIZE = CENTER_SIZE + 2 * HALO  # 384


def open_memmap(path):
    if not path.exists():
        raise FileNotFoundError(path)

    arr = tiff.memmap(str(path))

    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D volume, got shape={arr.shape}: {path}")

    return arr


def choose_center(shape):
    z, y, x = shape
    return z // 2, y // 2, x // 2


def crop_with_bounds(arr, center, size):
    zc, yc, xc = center
    half = size // 2

    z0 = max(0, zc - half)
    y0 = max(0, yc - half)
    x0 = max(0, xc - half)

    z1 = min(arr.shape[0], z0 + size)
    y1 = min(arr.shape[1], y0 + size)
    x1 = min(arr.shape[2], x0 + size)

    z0 = max(0, z1 - size)
    y0 = max(0, y1 - size)
    x0 = max(0, x1 - size)

    crop = np.asarray(arr[z0:z1, y0:y1, x0:x1])

    return crop, {
        "z0": int(z0), "z1": int(z1),
        "y0": int(y0), "y1": int(y1),
        "x0": int(x0), "x1": int(x1),
    }


def main():
    print("Opening gray memmap:")
    print(GRAY)
    gray = open_memmap(GRAY)

    print("Opening label memmap:")
    print(LABEL)
    label = open_memmap(LABEL)

    print("gray shape :", gray.shape, gray.dtype)
    print("label shape:", label.shape, label.dtype)

    if gray.shape != label.shape:
        raise RuntimeError(f"Shape mismatch: gray={gray.shape}, label={label.shape}")

    center = choose_center(gray.shape)

    print("chosen center:", center)
    print("processed tile size:", PROC_SIZE)
    print("center tile size:", CENTER_SIZE)
    print("halo:", HALO)

    gray_proc, bbox = crop_with_bounds(gray, center, PROC_SIZE)
    label_proc, bbox2 = crop_with_bounds(label, center, PROC_SIZE)

    if bbox != bbox2:
        raise RuntimeError("Bbox mismatch between gray and label.")

    start = HALO
    end = HALO + CENTER_SIZE

    gray_center = gray_proc[start:end, start:end, start:end]
    label_center = label_proc[start:end, start:end, start:end]

    print("gray_proc shape:", gray_proc.shape, gray_proc.dtype)
    print("label_proc shape:", label_proc.shape, label_proc.dtype)
    print("gray_center shape:", gray_center.shape, gray_center.dtype)
    print("label_center shape:", label_center.shape, label_center.dtype)

    print("label_proc nonzero voxels:", int(np.count_nonzero(label_proc)))
    print("label_center nonzero voxels:", int(np.count_nonzero(label_center)))
    print("label_proc unique labels approx:", int(len(np.unique(label_proc))) - 1)
    print("label_center unique labels approx:", int(len(np.unique(label_center))) - 1)

    tiff.imwrite(str(OUT / "gray_proc_384_with_halo.tif"), gray_proc.astype(np.uint16))
    tiff.imwrite(str(OUT / "label_proc_384_with_halo.tif"), label_proc.astype(np.uint16))

    tiff.imwrite(str(OUT / "gray_center_256.tif"), gray_center.astype(np.uint16))
    tiff.imwrite(str(OUT / "label_center_256.tif"), label_center.astype(np.uint16))

    meta = {
        "scan": "EFRGP01_01",
        "gray_source": str(GRAY),
        "label_source": str(LABEL),
        "full_shape": list(map(int, gray.shape)),
        "center_voxel": list(map(int, center)),
        "processed_tile_size": int(PROC_SIZE),
        "center_tile_size": int(CENTER_SIZE),
        "halo": int(HALO),
        "bbox_in_full_image": bbox,
        "gray_proc": str(OUT / "gray_proc_384_with_halo.tif"),
        "label_proc": str(OUT / "label_proc_384_with_halo.tif"),
        "gray_center": str(OUT / "gray_center_256.tif"),
        "label_center": str(OUT / "label_center_256.tif"),
        "label_proc_nonzero_voxels": int(np.count_nonzero(label_proc)),
        "label_center_nonzero_voxels": int(np.count_nonzero(label_center)),
        "label_proc_unique_labels": int(len(np.unique(label_proc))) - 1,
        "label_center_unique_labels": int(len(np.unique(label_center))) - 1,
    }

    with open(OUT / "tile_metadata.json", "w") as f:
        json.dump(meta, f, indent=2)

    print()
    print("Wrote pilot tile to:", OUT)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
