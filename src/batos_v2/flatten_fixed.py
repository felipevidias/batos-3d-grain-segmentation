"""Flatten-stage fix test (copy of src/batos/batos_preprocess.py; the frozen file is untouched).

Root cause (docs/CROP05_FORENSICS.md): flatten_volume_u8 sets bg_ref = median(blur[blur > 0]). When the openrec image
has large exactly-zero plateaux (exterior air + pores, 58-67 % of crop05), the float Gaussian tails over them are
counted as "positive", and bg_ref collapses to ~0.9 (186 on crop04). The selector then loses its offset.

This copy is identical to the frozen stage (same CLI, robust_u8, per-slice Gaussian, clip, robust p0.5-p99.5 rescale,
same output names) except for bg_ref:
  --variant A : bg_ref = median(blur[blur >= 1])                       (ignore sub-grey-level tails)
  --variant B : bg_ref = median(blur[Otsu foreground of the raw crop])  (e13's foreground, imported read-only)
  --variant frozen : bg_ref = median(blur[blur > 0])                    (reproduces the frozen stage; used as a check)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi

from .scale_granulo import load_e13


def read3d(path):
    arr = np.asarray(tiff.imread(str(path)))
    if arr.ndim == 2:
        arr = arr[np.newaxis, :, :]
    if arr.ndim == 4 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    if arr.ndim != 3:
        raise RuntimeError(f"Esperava volume 3D, mas shape={arr.shape}")
    return arr


def robust_u8(vol, p_low=0.5, p_high=99.5):
    vol = np.asarray(vol)
    if vol.dtype == np.uint8:
        return vol.copy()
    v = vol.astype(np.float32)
    nz = v[v > 0]
    if nz.size == 0:
        nz = v.ravel()
    lo = np.percentile(nz, p_low)
    hi = np.percentile(nz, p_high)
    if hi <= lo:
        hi = lo + 1.0
    out = (v - lo) * 255.0 / (hi - lo)
    return np.clip(out, 0, 255).astype(np.uint8)


def flatten_volume_u8(vol_u8, lambda_xy, sigma_z, variant, fg=None):
    vol = vol_u8.astype(np.float32)
    if sigma_z > 0:
        bg = ndi.gaussian_filter(vol, sigma=(sigma_z, lambda_xy, lambda_xy))
    else:
        bg = np.empty_like(vol)
        for z in range(vol.shape[0]):
            bg[z] = ndi.gaussian_filter(vol[z], sigma=lambda_xy)
    if variant == "frozen":
        bg_ref = np.median(bg[bg > 0]) if np.any(bg > 0) else np.median(bg)
    elif variant == "A":
        bg_ref = np.median(bg[bg >= 1]) if np.any(bg >= 1) else np.median(bg)
    elif variant == "B":
        if fg is None or fg.shape != bg.shape or not fg.any():
            raise RuntimeError("variant B needs a non-empty foreground of the same shape")
        bg_ref = np.median(bg[fg])
    else:
        raise ValueError(variant)
    flat = vol - bg + bg_ref
    flat = np.clip(flat, 0, None)
    return robust_u8(flat, p_low=0.5, p_high=99.5), float(bg_ref)


def main():
    ap = argparse.ArgumentParser(description="flatten stage with a specimen-only offset (fix test)")
    ap.add_argument("input_tif")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default=None)
    ap.add_argument("--lambda-xy", type=float, default=15.0)
    ap.add_argument("--sigma-z", type=float, default=0.0)
    ap.add_argument("--variant", choices=("A", "B", "frozen"), required=True)
    ap.add_argument("--raw-gray", help="raw crop (variant B: e13 Otsu foreground)")
    args = ap.parse_args()

    inp = Path(args.input_tif)
    prefix = args.prefix or inp.stem
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    vol = read3d(inp)
    orig8 = robust_u8(vol)
    fg = None
    if args.variant == "B":
        e13 = load_e13()
        fg, _ = e13.make_foreground(e13.robust_u8(np.squeeze(read3d(args.raw_gray))))
    flat, bg_ref = flatten_volume_u8(orig8, args.lambda_xy, args.sigma_z, args.variant, fg)

    lambda_tag = str(args.lambda_xy).replace(".", "p")
    orig_path = out / f"{prefix}_original_8u.tif"
    flat_path = out / f"{prefix}_flatten_lambda{lambda_tag}_8u.tif"
    tiff.imwrite(str(orig_path), orig8)
    tiff.imwrite(str(flat_path), flat)
    (out / "flatten_fixed_info.json").write_text(json.dumps(
        {"variant": args.variant, "bg_ref": bg_ref, "raw_gray": args.raw_gray, "input": str(inp)}, indent=1),
        encoding="utf-8")
    print(f"[flatten-fixed] variant={args.variant} bg_ref={bg_ref:.3f} flat mean {flat.mean():.3f}")


if __name__ == "__main__":
    main()
