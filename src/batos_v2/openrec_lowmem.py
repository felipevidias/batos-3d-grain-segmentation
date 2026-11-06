"""Low-memory, bit-exact re-implementation of the frozen openrec v2 stage.

Frozen script: linux_import/zsh_recovered/code__segmentation__batos__batos_preprocess_opening_reconstruction_v2.py
(unchanged; this module only re-implements it). The frozen erosion `skimage.morphology.erosion(img, ball(r))` calls
scipy grey_erosion with a 3-D footprint; scipy allocates an offset table of about nonzero(ball) x (2r+1)^3 entries,
independent of the image size (~3.3 GB at r = 15, ~18 GB at r = 20, MemoryError at r >= 27 on 32 GB).

Here the same erosion is computed exactly as a minimum over z-shifts of 2-D erosions by the z-slices of ball(r):
  erosion_ball(f)(z, y, x) = min_{dz} erosion_{disk_dz}(f(z + dz, ., .))(y, x),
with scipy's 'reflect' boundary handled per axis (the z axis by symmetric padding, which equals scipy 'reflect';
y/x by grey_erosion itself). The min over a union of footprints is exact, so the result is bit-identical; offset
memory drops to the 2-D tables. Everything else (robust_u8, reconstruction by dilation, flattening, rescaling, output
names and report) is copied from the frozen script verbatim. CLI identical to the frozen script.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scipy.ndimage as ndi
import tifffile as tiff
from skimage.morphology import ball, reconstruction


def read3d(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise RuntimeError(f"Esperado volume 3D, veio {arr.shape}: {path}")
    return arr


def robust_u8(vol):
    v = vol.astype(np.float32)
    p1, p99 = np.percentile(v, [1, 99])
    v = (v - p1) / (p99 - p1 + 1e-8)
    v = np.clip(v, 0, 1)
    return (255 * v).astype(np.uint8), float(p1), float(p99)


rescale_u8 = robust_u8


def erosion_ball_lowmem(img: np.ndarray, radius: int) -> np.ndarray:
    se = ball(radius).astype(bool)
    r = radius
    nz = img.shape[0]
    padded = np.pad(img, ((r, r), (0, 0), (0, 0)), mode="symmetric")  # == scipy 'reflect' along z
    out = np.full(img.shape, np.iinfo(img.dtype).max, dtype=img.dtype)
    for d in range(0, r + 1):
        fp = se[r + d]
        if not fp.any():
            continue
        e2 = ndi.grey_erosion(padded, footprint=fp[None, :, :], mode="reflect")
        np.minimum(out, e2[r + d: r + d + nz], out=out)
        if d:
            np.minimum(out, e2[r - d: r - d + nz], out=out)
        del e2
    return out


def opening_by_reconstruction(img_u8, radius):
    marker = erosion_ball_lowmem(img_u8, radius)
    rec = reconstruction(marker, img_u8, method="dilation")
    return rec.astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_tif")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--radius", type=int, required=True)
    args = ap.parse_args()

    inp = Path(args.input_tif)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[openrec-v2 lowmem] input :", inp)
    print("[openrec-v2 lowmem] radius:", args.radius)

    vol = read3d(inp)
    img_u8, p1, p99 = robust_u8(vol)
    opened = opening_by_reconstruction(img_u8, args.radius)
    flattened = img_u8.astype(np.float32) - opened + float(np.median(opened))
    flattened_u8, f1, f99 = rescale_u8(flattened)
    tophat = img_u8.astype(np.float32) - opened
    tophat_u8, t1, t99 = rescale_u8(tophat)

    tag = f"r{args.radius}"
    out_original = out_dir / f"{args.prefix}_original_8u.tif"
    out_bg = out_dir / f"{args.prefix}_debug_opening_reconstruction_background_{tag}_8u.tif"
    out_flat = out_dir / f"{args.prefix}_flatten_openrec_{tag}_8u.tif"
    out_tophat = out_dir / f"{args.prefix}_openrec_tophat_{tag}_8u.tif"
    out_report = out_dir / f"{args.prefix}_openrec_v2_report.json"
    tiff.imwrite(str(out_original), img_u8)
    tiff.imwrite(str(out_bg), np.clip(opened, 0, 255).astype(np.uint8))
    tiff.imwrite(str(out_flat), flattened_u8)
    tiff.imwrite(str(out_tophat), tophat_u8)
    report = {"input": str(inp), "shape": list(vol.shape), "dtype": str(vol.dtype), "radius": int(args.radius),
              "input_p1": p1, "input_p99": p99, "flattened_p1": f1, "flattened_p99": f99, "tophat_p1": t1,
              "tophat_p99": t99, "original_8u": str(out_original), "background_debug": str(out_bg),
              "flattened_openrec": str(out_flat), "openrec_tophat": str(out_tophat),
              "implementation": "batos_v2.openrec_lowmem (z-decomposed exact ball erosion)"}
    out_report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("wrote:", out_flat)


if __name__ == "__main__":
    main()
