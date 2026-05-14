#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.feature import peak_local_max
from skimage.segmentation import watershed


def read3d(path):
    arr = np.asarray(tiff.imread(str(path)))
    if arr.ndim == 4 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D volume, got shape={arr.shape} for {path}")
    return arr


def robust_u8(vol):
    if vol.dtype == np.uint8:
        return vol.copy()

    v = vol.astype(np.float32)
    nz = v[v > 0]
    if nz.size == 0:
        nz = v.ravel()

    lo = np.percentile(nz, 0.5)
    hi = np.percentile(nz, 99.5)
    if hi <= lo:
        hi = lo + 1.0

    return np.clip((v - lo) * 255.0 / (hi - lo), 0, 255).astype(np.uint8)


def make_ball(radius):
    r = int(radius)
    if r <= 0:
        return np.ones((1, 1, 1), dtype=bool)
    zz, yy, xx = np.ogrid[-r:r+1, -r:r+1, -r:r+1]
    return (xx * xx + yy * yy + zz * zz) <= r * r


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


def overlay_labels(gray_u8, labels, alpha=0.55):
    rgb = np.repeat(gray_u8[..., None], 3, axis=-1).astype(np.float32)
    col = colorize_labels(labels).astype(np.float32)
    mask = labels > 0

    out = rgb.copy()
    out[mask] = (1.0 - alpha) * rgb[mask] + alpha * col[mask]
    return np.clip(out, 0, 255).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--image", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="otsu_ws")
    ap.add_argument("--threshold-offset", type=float, default=-10)
    ap.add_argument("--min-distance", type=int, default=8)
    ap.add_argument("--close-radius", type=int, default=1)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gray = robust_u8(read3d(args.image))

    otsu = threshold_otsu(gray)
    thr = int(np.clip(otsu + args.threshold_offset, 0, 255))

    mask = gray >= thr

    if args.close_radius > 0:
        mask = ndi.binary_closing(mask, structure=make_ball(args.close_radius))

    dist = ndi.distance_transform_edt(mask)

    coords = peak_local_max(
        dist,
        min_distance=args.min_distance,
        labels=mask,
        exclude_border=False
    )

    markers = np.zeros_like(gray, dtype=np.int32)
    for i, (z, y, x) in enumerate(coords, start=1):
        markers[z, y, x] = i

    labels = watershed(
        -dist,
        markers=markers,
        mask=mask,
        connectivity=1
    )

    ids = np.unique(labels)
    ids = ids[ids > 0]

    labels_out = out_dir / f"{args.prefix}_labels.tif"
    markers_out = out_dir / f"{args.prefix}_markers.tif"
    mask_out = out_dir / f"{args.prefix}_mask.tif"
    colored_out = out_dir / f"{args.prefix}_colored_labels_rgb.tif"
    overlay_out = out_dir / f"{args.prefix}_overlay_rgb.tif"

    dtype = np.uint16 if int(labels.max()) <= 65535 else np.uint32

    tiff.imwrite(labels_out, labels.astype(dtype))
    tiff.imwrite(markers_out, markers.astype(np.uint16 if markers.max() <= 65535 else np.uint32))
    tiff.imwrite(mask_out, (mask.astype(np.uint8) * 255))
    tiff.imwrite(colored_out, colorize_labels(labels))
    tiff.imwrite(overlay_out, overlay_labels(gray, labels))

    print("[otsu-ws] image shape     :", gray.shape)
    print("[otsu-ws] otsu threshold  :", otsu)
    print("[otsu-ws] used threshold  :", thr)
    print("[otsu-ws] foreground voxels:", int(mask.sum()))
    print("[otsu-ws] markers         :", len(coords))
    print("[otsu-ws] labels          :", len(ids))
    print("[otsu-ws] wrote labels    :", labels_out)
    print("[otsu-ws] wrote overlay   :", overlay_out)
    print("[otsu-ws] open:")
    print(f"gmic {overlay_out} a z")
    print(f"gmic {colored_out} a z")


if __name__ == "__main__":
    main()
