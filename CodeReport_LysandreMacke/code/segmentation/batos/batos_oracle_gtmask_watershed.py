#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi
from skimage.segmentation import watershed


def read3d(path):
    arr = np.asarray(tiff.imread(str(path)))
    if arr.ndim == 4 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D volume, got {arr.shape}")
    return arr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gt", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="oracle_gtmask")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gt = read3d(args.gt)
    markers = read3d(args.markers)

    mask = gt > 0
    marker_mask = markers > 0

    marker_lab, n_markers = ndi.label(marker_mask)

    # Remove markers outside GT mask
    marker_lab[~mask] = 0
    marker_lab, n_markers_inside = ndi.label(marker_lab > 0)

    print("[oracle] GT foreground voxels:", int(mask.sum()))
    print("[oracle] connected markers inside GT:", int(n_markers_inside))

    dist = ndi.distance_transform_edt(mask)

    labels = watershed(
        -dist,
        markers=marker_lab.astype(np.int32),
        mask=mask,
        connectivity=1
    )

    out = out_dir / f"{args.prefix}_labels.tif"
    tiff.imwrite(out, labels.astype(np.uint16))

    ids = np.unique(labels)
    ids = ids[ids > 0]

    print("[oracle] output labels:", len(ids))
    print("[oracle] wrote:", out)


if __name__ == "__main__":
    main()
