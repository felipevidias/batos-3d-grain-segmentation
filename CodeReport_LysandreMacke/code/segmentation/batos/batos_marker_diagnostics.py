#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi
import csv


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
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gt = read3d(args.gt)
    markers = read3d(args.markers)

    marker_mask = markers > 0
    marker_lab, n_markers = ndi.label(marker_mask)

    ids = np.arange(1, n_markers + 1)
    rows = []

    markers_inside_gt = 0
    markers_on_bg = 0
    marker_voxels_on_gt = int(np.logical_and(marker_mask, gt > 0).sum())
    marker_voxels_on_bg = int(np.logical_and(marker_mask, gt == 0).sum())

    gt_hit = set()

    for mid in ids:
        m = marker_lab == mid
        touched = gt[m]
        touched = touched[touched > 0]

        if touched.size == 0:
            markers_on_bg += 1
            rows.append([mid, int(m.sum()), 0, 0, "background"])
            continue

        vals, counts = np.unique(touched, return_counts=True)
        j = int(np.argmax(counts))
        dominant_gt = int(vals[j])
        overlap = int(counts[j])

        markers_inside_gt += 1
        gt_hit.add(dominant_gt)
        rows.append([mid, int(m.sum()), dominant_gt, overlap, "inside_gt"])

    with open(out_dir / "marker_diagnostics.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["marker_id", "marker_voxels", "dominant_gt_id", "overlap_voxels", "status"])
        w.writerows(rows)

    with open(out_dir / "marker_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["marker_components_connected", n_markers])
        w.writerow(["markers_inside_gt", markers_inside_gt])
        w.writerow(["markers_on_background", markers_on_bg])
        w.writerow(["marker_voxels_on_gt", marker_voxels_on_gt])
        w.writerow(["marker_voxels_on_background", marker_voxels_on_bg])
        w.writerow(["unique_gt_objects_hit_by_markers", len(gt_hit)])

    print("marker_components_connected:", n_markers)
    print("markers_inside_gt:", markers_inside_gt)
    print("markers_on_background:", markers_on_bg)
    print("marker_voxels_on_gt:", marker_voxels_on_gt)
    print("marker_voxels_on_background:", marker_voxels_on_bg)
    print("unique_gt_objects_hit_by_markers:", len(gt_hit))
    print("wrote:", out_dir)


if __name__ == "__main__":
    main()
