#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.segmentation import watershed


def read3d(path):
    arr = np.asarray(tiff.imread(str(path)))
    if arr.ndim == 4 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D volume, got {arr.shape} for {path}")
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
        hi = lo + 1
    out = (v - lo) * 255.0 / (hi - lo)
    return np.clip(out, 0, 255).astype(np.uint8)


def make_ball(radius):
    r = int(radius)
    if r <= 0:
        return np.ones((1, 1, 1), dtype=bool)
    zz, yy, xx = np.ogrid[-r:r+1, -r:r+1, -r:r+1]
    return (xx*xx + yy*yy + zz*zz) <= r*r


def colorize_labels(labels):
    labels = labels.astype(np.int64)
    out = np.zeros(labels.shape + (3,), dtype=np.uint8)

    ids = np.unique(labels)
    ids = ids[ids > 0]

    for lab in ids:
        # hash determinístico simples
        r = (lab * 37 + 53) % 255
        g = (lab * 73 + 91) % 255
        b = (lab * 109 + 17) % 255
        out[labels == lab] = (r, g, b)

    return out


def overlay_labels(gray_u8, labels, alpha=0.55):
    rgb_gray = np.repeat(gray_u8[..., None], 3, axis=-1).astype(np.float32)
    col = colorize_labels(labels).astype(np.float32)
    mask = labels > 0

    out = rgb_gray.copy()
    out[mask] = (1.0 - alpha) * rgb_gray[mask] + alpha * col[mask]
    return np.clip(out, 0, 255).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--image", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="batos_distance")
    ap.add_argument("--threshold-offset", type=float, default=-10.0)
    ap.add_argument("--close-radius", type=int, default=2)
    ap.add_argument("--fill-holes", action="store_true")
    ap.add_argument("--force-markers-foreground", action="store_true")
    ap.add_argument("--force-radius", type=int, default=1)
    ap.add_argument("--min-object-size", type=int, default=0)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    img = read3d(args.image)
    gray = robust_u8(img)

    marker_raw = read3d(args.markers)
    marker_mask = marker_raw > 0
    marker_labels, n_markers = ndi.label(marker_mask)

    otsu = threshold_otsu(gray)
    thr = int(np.clip(otsu + args.threshold_offset, 0, 255))

    # Grãos são foreground claro no crop.
    fg = gray >= thr

    if args.close_radius > 0:
        ball = make_ball(args.close_radius)
        fg = ndi.binary_closing(fg, structure=ball)

    if args.fill_holes:
        # fill holes 3D; em alguns volumes pode ser agressivo.
        fg = ndi.binary_fill_holes(fg)

    if args.min_object_size > 0:
        lab, n = ndi.label(fg)
        counts = np.bincount(lab.ravel())
        keep = np.zeros_like(counts, dtype=bool)
        keep[0] = False
        keep[counts >= args.min_object_size] = True
        fg = keep[lab]

    if args.force_markers_foreground:
        if args.force_radius > 0:
            marker_fg = ndi.binary_dilation(marker_mask, structure=make_ball(args.force_radius))
        else:
            marker_fg = marker_mask
        fg = np.logical_or(fg, marker_fg)

    # Remove markers fora do foreground final.
    marker_labels[~fg] = 0
    marker_labels, n_markers_inside = ndi.label(marker_labels > 0)

    dist = ndi.distance_transform_edt(fg)

    labels = watershed(
        -dist,
        markers=marker_labels.astype(np.int32),
        mask=fg,
        connectivity=1
    )

    ids = np.unique(labels)
    ids = ids[ids > 0]

    labels_out = out_dir / f"{args.prefix}_labels.tif"
    mask_out = out_dir / f"{args.prefix}_foreground_mask.tif"
    markers_out = out_dir / f"{args.prefix}_markers_inside_fg.tif"
    colored_out = out_dir / f"{args.prefix}_colored_labels_rgb.tif"
    overlay_out = out_dir / f"{args.prefix}_overlay_rgb.tif"

    dtype = np.uint16 if int(labels.max()) <= 65535 else np.uint32
    tiff.imwrite(labels_out, labels.astype(dtype))
    tiff.imwrite(mask_out, (fg.astype(np.uint8) * 255))
    tiff.imwrite(markers_out, marker_labels.astype(np.uint16))
    tiff.imwrite(colored_out, colorize_labels(labels))
    tiff.imwrite(overlay_out, overlay_labels(gray, labels))

    print("[batos-distance] image shape          :", gray.shape)
    print("[batos-distance] otsu threshold       :", otsu)
    print("[batos-distance] used threshold       :", thr)
    print("[batos-distance] foreground voxels    :", int(fg.sum()))
    print("[batos-distance] raw marker components:", int(n_markers))
    print("[batos-distance] markers inside fg    :", int(n_markers_inside))
    print("[batos-distance] final labels         :", int(len(ids)))
    print("[batos-distance] wrote labels         :", labels_out)
    print("[batos-distance] wrote overlay        :", overlay_out)
    print("[batos-distance] open:")
    print(f"gmic {overlay_out} a z")
    print(f"gmic {colored_out} a z")


if __name__ == "__main__":
    main()
