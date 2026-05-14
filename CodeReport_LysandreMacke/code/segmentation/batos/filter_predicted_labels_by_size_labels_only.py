#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi

def read3d(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise RuntimeError(f"expected 3D image, got shape={arr.shape}")
    return arr

def relabel_sequential(labels):
    ids = np.unique(labels)
    ids = ids[ids > 0]
    out = np.zeros(labels.shape, dtype=np.uint16)
    for new_id, old_id in enumerate(ids, start=1):
        out[labels == old_id] = new_id
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="filtered")
    ap.add_argument("--min-size", type=int, required=True)
    ap.add_argument("--relabel", action="store_true")
    args = ap.parse_args()

    labels = read3d(args.labels)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ids, counts = np.unique(labels, return_counts=True)
    keep_ids = ids[(ids > 0) & (counts >= args.min_size)]

    keep = np.zeros(int(ids.max()) + 1 if ids.size else 1, dtype=bool)
    keep[keep_ids] = True

    filtered = np.where((labels > 0) & keep[labels], labels, 0).astype(np.uint16)

    before = int(len(ids[ids > 0]))
    if args.relabel:
        filtered = relabel_sequential(filtered)

    after = int(len(np.unique(filtered)) - 1)

    labels_out = out_dir / f"{args.prefix}_min{args.min_size}_labels.tif"
    tiff.imwrite(labels_out, filtered)

    print("[filter-size-labels-only] input labels :", args.labels)
    print("[filter-size-labels-only] min-size     :", args.min_size)
    print("[filter-size-labels-only] before labels:", before)
    print("[filter-size-labels-only] after labels :", after)
    print("[filter-size-labels-only] removed      :", before - after)
    print("[filter-size-labels-only] wrote labels :", labels_out)

if __name__ == "__main__":
    main()
