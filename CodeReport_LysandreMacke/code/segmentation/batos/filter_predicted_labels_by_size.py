#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import tifffile as tiff


def read3d(path):
    arr = np.asarray(tiff.imread(str(path)))
    if arr.ndim == 4 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D label volume, got {arr.shape}")
    return arr


def relabel_sequential(labels):
    ids = np.unique(labels)
    ids = ids[ids > 0]

    out = np.zeros_like(labels, dtype=np.uint32)
    for new_id, old_id in enumerate(ids, start=1):
        out[labels == old_id] = new_id

    return out


def colorize_labels(labels):
    labels = labels.astype(np.int64)
    out = np.zeros(labels.shape + (3,), dtype=np.uint8)

    ids = np.unique(labels)
    ids = ids[ids > 0]

    for lab in ids:
        r = (lab * 37 + 53) % 255
        g = (lab * 73 + 91) % 255
        b = (lab * 109 + 17) % 255
        out[labels == lab] = (r, g, b)

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="filtered")
    ap.add_argument("--min-size", type=int, required=True)
    ap.add_argument("--relabel", action="store_true")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    labels = read3d(args.labels).astype(np.int64)

    ids, counts = np.unique(labels, return_counts=True)
    keep_ids = ids[(ids > 0) & (counts >= args.min_size)]

    keep = np.isin(labels, keep_ids)
    out = labels.copy()
    out[~keep] = 0

    before = len(ids[ids > 0])
    after = len(keep_ids)

    if args.relabel:
        out = relabel_sequential(out)

    dtype = np.uint16 if int(out.max()) <= 65535 else np.uint32

    labels_out = out_dir / f"{args.prefix}_min{args.min_size}_labels.tif"
    color_out = out_dir / f"{args.prefix}_min{args.min_size}_colored_labels_rgb.tif"

    tiff.imwrite(labels_out, out.astype(dtype))
    tiff.imwrite(color_out, colorize_labels(out))

    print("[filter-size] input labels :", args.labels)
    print("[filter-size] min-size     :", args.min_size)
    print("[filter-size] before labels:", before)
    print("[filter-size] after labels :", after)
    print("[filter-size] removed      :", before - after)
    print("[filter-size] wrote labels :", labels_out)
    print("[filter-size] wrote color  :", color_out)


if __name__ == "__main__":
    main()
