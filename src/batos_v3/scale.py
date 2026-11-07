"""BA-TOS v3 scale source: r_granulo = argmax of the binary pattern spectrum of the specimen-only Otsu foreground.

GT-free, no free parameter beyond the radius range r = 2..40 fixed by the design (docs/V3_SCALE_FREE.md):

  pass 1  u8_0  = robust p1-p99 rescale of the crop;  fg0 = fill_holes(u8_0 > Otsu(u8_0)) (frozen otsu_u8 convention: dark <= t)
          r0    = argmax_r PS_fg0(r)                                (first-pass estimate)
  mask    specimen = fill_holes(largest 6-connected component of closing(fg0, ball(2 r0)))
  pass 2  u8_1  = robust p1-p99 rescale computed on specimen voxels; fg1 = fill_holes((u8_1 > Otsu(u8_1[specimen])) & specimen)
          r_granulo = argmax_r PS_fg1(r)
  openrec radius = max(8, round(1.5 r_granulo))

Pattern spectrum: V(r) = volume of the opening of fg by the ball |d|^2 <= r^2 (skimage ball), computed exactly with
Euclidean distance transforms (erosion = {EDT_fg > r}; dilation of that set = {EDT to it <= r}); PS(r) = V(r) - V(r+1),
r = 2..40. Voxels outside the crop are not background (grains cut by a face are not shrunk by the face); voxels
outside the specimen are background (pass 2). fill_holes (parameter-free) removes enclosed noise cavities inside grains,
which otherwise cap the EDT and move the peak to r = 7 on every dev crop (design iteration 1, EXPERIMENT_LOG). Closing: dilation {EDT to fg0 <= R}, then erosion {EDT to the
complement > R}, again with the crop outside not treated as background.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scipy.ndimage as ndi
import tifffile

R_MIN, R_MAX = 2, 40
CONN6 = ndi.generate_binary_structure(3, 1)


def robust_u8(values: np.ndarray, ref: np.ndarray) -> np.ndarray:
    p1, p99 = np.percentile(ref.astype(np.float32), [1, 99])
    v = (values.astype(np.float32) - p1) / (p99 - p1 + 1e-8)
    return (255 * np.clip(v, 0, 1)).astype(np.uint8)


def otsu_u8(values: np.ndarray) -> int:
    """Otsu on a uint8 histogram (same formula as the frozen batos_background.otsu_u8; returns t, fg = v >= t+1)."""
    hist = np.bincount(values.ravel(), minlength=256).astype(np.float64)
    p = hist / max(1.0, hist.sum())
    omega = np.cumsum(p)
    mu = np.cumsum(p * np.arange(256))
    denom = omega * (1.0 - omega)
    denom[denom <= 1e-12] = np.nan
    sigma = (mu[-1] * omega - mu) ** 2 / denom
    return int(np.nanargmax(sigma))


def pattern_spectrum(fg: np.ndarray) -> tuple[int, dict, dict]:
    edt = ndi.distance_transform_edt(fg)
    vols = {}
    for r in range(R_MIN, R_MAX + 2):
        core = edt > r
        vols[r] = int((ndi.distance_transform_edt(~core) <= r).sum()) if core.any() else 0
    ps = {r: vols[r] - vols[r + 1] for r in range(R_MIN, R_MAX + 1)}
    r_peak = max(ps, key=lambda r: (ps[r], -r))  # argmax; ties -> smallest radius
    return int(r_peak), ps, vols


def closing(mask: np.ndarray, radius: float) -> np.ndarray:
    if mask.all():
        return mask.copy()
    dil = ndi.distance_transform_edt(~mask) <= radius
    if dil.all():
        return dil
    return ndi.distance_transform_edt(dil) > radius


def specimen_mask(fg0: np.ndarray, r0: int) -> np.ndarray:
    closed = closing(fg0, 2 * r0)
    lab, n = ndi.label(closed, structure=CONN6)
    if n == 0:
        raise RuntimeError("empty specimen mask")
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    return ndi.binary_fill_holes(lab == int(np.argmax(sizes)))


def estimate(gray: np.ndarray) -> tuple[dict, np.ndarray]:
    gray = np.squeeze(gray)
    u0 = robust_u8(gray, gray)
    t0 = otsu_u8(u0)
    fg0 = ndi.binary_fill_holes(u0 > t0)  # enclosed noise cavities would cap the EDT (iteration 1)
    r0, ps0, _ = pattern_spectrum(fg0)
    spec = specimen_mask(fg0, r0)
    u1 = robust_u8(gray, gray[spec])
    t1 = otsu_u8(u1[spec])
    fg1 = ndi.binary_fill_holes((u1 > t1) & spec) & spec
    r, ps1, vols1 = pattern_spectrum(fg1)
    res = {"r0": r0, "otsu_u8_pass1": t0, "fg_fraction_pass1": float(fg0.mean()),
           "specimen_fraction": float(spec.mean()), "otsu_u8_specimen": t1,
           "fg_fraction_in_specimen": float(fg1.sum() / max(1, spec.sum())),
           "r_granulo": r, "openrec_radius": int(max(8, round(1.5 * r))),
           "radius_range": [R_MIN, R_MAX], "gt_used": False,
           "pattern_spectrum_pass1": ps0, "pattern_spectrum": ps1, "opened_volume": vols1}
    return res, spec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="BA-TOS v3 scale: r_granulo and specimen mask (GT-free)")
    ap.add_argument("gray", type=Path)
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--specimen-tif", type=Path, required=True)
    args = ap.parse_args(argv)
    res, spec = estimate(tifffile.imread(args.gray))
    res["gray"] = str(args.gray)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(res, indent=1), encoding="utf-8")
    tifffile.imwrite(args.specimen_tif, spec.astype(np.uint8), compression="zlib")
    print(json.dumps({k: v for k, v in res.items() if not isinstance(v, dict)}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
