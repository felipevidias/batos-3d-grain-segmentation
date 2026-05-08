#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi

def read3d(path):
    arr = tiff.imread(str(path))
    arr = np.asarray(arr)
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

def flatten_volume_u8(vol_u8, lambda_xy=15.0, sigma_z=0.0):
    vol = vol_u8.astype(np.float32)
    if sigma_z > 0:
        bg = ndi.gaussian_filter(vol, sigma=(sigma_z, lambda_xy, lambda_xy))
    else:
        bg = np.empty_like(vol)
        for z in range(vol.shape[0]):
            bg[z] = ndi.gaussian_filter(vol[z], sigma=lambda_xy)
    bg_ref = np.median(bg[bg > 0]) if np.any(bg > 0) else np.median(bg)
    flat = vol - bg + bg_ref
    flat = np.clip(flat, 0, None)
    return robust_u8(flat, p_low=0.5, p_high=99.5)

def main():
    ap = argparse.ArgumentParser(description="v11: flatten conservador lambda=15, sem core_enhanced.")
    ap.add_argument("input_tif")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default=None)
    ap.add_argument("--lambda-xy", type=float, default=15.0)
    ap.add_argument("--sigma-z", type=float, default=0.0)
    args = ap.parse_args()

    inp = Path(args.input_tif)
    prefix = args.prefix or inp.stem
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    vol = read3d(inp)
    orig8 = robust_u8(vol)
    flat = flatten_volume_u8(orig8, lambda_xy=args.lambda_xy, sigma_z=args.sigma_z)

    lambda_tag = str(args.lambda_xy).replace(".", "p")
    orig_path = out / f"{prefix}_original_8u.tif"
    flat_path = out / f"{prefix}_flatten_lambda{lambda_tag}_8u.tif"

    tiff.imwrite(str(orig_path), orig8)
    tiff.imwrite(str(flat_path), flat)

    print(f"[flatten-v11] input       = {inp}")
    print(f"[flatten-v11] out_dir     = {out}")
    print(f"[flatten-v11] lambda_xy   = {args.lambda_xy}")
    print(f"[flatten-v11] original    = {orig_path}")
    print(f"[flatten-v11] flattened   = {flat_path}")
    print(f"[flatten-v11] flat stats  = min {flat.min()}, max {flat.max()}, mean {flat.mean():.3f}, std {flat.std():.3f}")

if __name__ == "__main__":
    main()
