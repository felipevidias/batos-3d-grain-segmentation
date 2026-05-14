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


def inspect_tiff(path):
    if not path.exists():
        raise FileNotFoundError(path)

    with tiff.TiffFile(str(path)) as tif:
        series = tif.series[0]
        shape = tuple(series.shape)
        dtype = series.dtype
        axes = series.axes

    if len(shape) != 3:
        raise RuntimeError(f"Expected 3D TIFF, got shape={shape}: {path}")

    return shape, dtype, axes


def choose_center(shape):
    z, y, x = shape
    return z // 2, y // 2, x // 2


def bbox_from_center(shape, center, size):
    zc, yc, xc = center
    half = size // 2

    z0 = max(0, zc - half)
    y0 = max(0, yc - half)
    x0 = max(0, xc - half)

    z1 = min(shape[0], z0 + size)
    y1 = min(shape[1], y0 + size)
    x1 = min(shape[2], x0 + size)

    z0 = max(0, z1 - size)
    y0 = max(0, y1 - size)
    x0 = max(0, x1 - size)

    return {
        "z0": int(z0), "z1": int(z1),
        "y0": int(y0), "y1": int(y1),
        "x0": int(x0), "x1": int(x1),
    }


def read_tile_by_pages(path, bbox):
    """
    Lê somente z0:z1 e recorta y/x por página.
    Não usa memmap. Funciona em partição somente leitura.
    """
    z0, z1 = bbox["z0"], bbox["z1"]
    y0, y1 = bbox["y0"], bbox["y1"]
    x0, x1 = bbox["x0"], bbox["x1"]

    slices = []

    with tiff.TiffFile(str(path)) as tif:
        pages = tif.pages

        if len(pages) < z1:
            raise RuntimeError(
                f"TIFF tem {len(pages)} páginas, mas preciso até z={z1}: {path}"
            )

        for z in range(z0, z1):
            if (z - z0) % 25 == 0:
                print(f"  reading z {z}/{z1} from {path.name}")

            page = pages[z].asarray()
            crop2d = page[y0:y1, x0:x1]
            slices.append(crop2d)

    arr = np.stack(slices, axis=0)
    return arr


def main():
    print("Inspecting gray:")
    gray_shape, gray_dtype, gray_axes = inspect_tiff(GRAY)
    print("gray path :", GRAY)
    print("gray shape:", gray_shape)
    print("gray dtype:", gray_dtype)
    print("gray axes :", gray_axes)

    print()
    print("Inspecting label:")
    label_shape, label_dtype, label_axes = inspect_tiff(LABEL)
    print("label path :", LABEL)
    print("label shape:", label_shape)
    print("label dtype:", label_dtype)
    print("label axes :", label_axes)

    if gray_shape != label_shape:
        raise RuntimeError(f"Shape mismatch: gray={gray_shape}, label={label_shape}")

    center = choose_center(gray_shape)
    bbox = bbox_from_center(gray_shape, center, PROC_SIZE)

    print()
    print("chosen center:", center)
    print("processed tile size:", PROC_SIZE)
    print("center tile size:", CENTER_SIZE)
    print("halo:", HALO)
    print("bbox:", bbox)

    print()
    print("Reading gray tile safely...")
    gray_proc = read_tile_by_pages(GRAY, bbox)

    print()
    print("Reading label tile safely...")
    label_proc = read_tile_by_pages(LABEL, bbox)

    start = HALO
    end = HALO + CENTER_SIZE

    gray_center = gray_proc[start:end, start:end, start:end]
    label_center = label_proc[start:end, start:end, start:end]

    print()
    print("gray_proc shape:", gray_proc.shape, gray_proc.dtype)
    print("label_proc shape:", label_proc.shape, label_proc.dtype)
    print("gray_center shape:", gray_center.shape, gray_center.dtype)
    print("label_center shape:", label_center.shape, label_center.dtype)

    label_proc_unique = int(len(np.unique(label_proc)) - 1)
    label_center_unique = int(len(np.unique(label_center)) - 1)

    print("label_proc nonzero voxels:", int(np.count_nonzero(label_proc)))
    print("label_center nonzero voxels:", int(np.count_nonzero(label_center)))
    print("label_proc unique labels:", label_proc_unique)
    print("label_center unique labels:", label_center_unique)

    tiff.imwrite(str(OUT / "gray_proc_384_with_halo.tif"), gray_proc.astype(np.uint16))
    tiff.imwrite(str(OUT / "label_proc_384_with_halo.tif"), label_proc.astype(np.uint16))

    tiff.imwrite(str(OUT / "gray_center_256.tif"), gray_center.astype(np.uint16))
    tiff.imwrite(str(OUT / "label_center_256.tif"), label_center.astype(np.uint16))

    meta = {
        "scan": "EFRGP01_01",
        "gray_source": str(GRAY),
        "label_source": str(LABEL),
        "full_shape": list(map(int, gray_shape)),
        "full_dtype_gray": str(gray_dtype),
        "full_dtype_label": str(label_dtype),
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
        "label_proc_unique_labels": label_proc_unique,
        "label_center_unique_labels": label_center_unique,
        "read_mode": "page_by_page_no_memmap",
    }

    with open(OUT / "tile_metadata.json", "w") as f:
        json.dump(meta, f, indent=2)

    print()
    print("Wrote pilot tile to:", OUT)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
