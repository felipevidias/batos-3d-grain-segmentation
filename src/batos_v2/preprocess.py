"""Stage 1: GT-free scale estimate, opening-by-reconstruction flattening, lambda-1 selector image.

Recovered:
- scale estimate: exact port of linux_patches/.../e13_estimate_scale.py (BG1)
- lambda flattening: exact port of src/batos/batos_preprocess.py flatten_volume_u8 (BG6)
Reconstructed (the original batos_preprocess_opening_reconstruction_v2.py is lost):
- BG2 radius, BG3 footprint, BG4 background use, BG5 polarity, BG7 input dtype

Outputs (E13 naming):
  <prefix>_flatten_openrec_r<R>_8u.tif   watershed / control image ("FLAT_IMAGE")
  <prefix>_flatten_lambda1p0_8u.tif      selector image for ToS + boundary step
  <prefix>_preprocess.json               parameters and statistics
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage as ndi
from skimage import morphology
from skimage.filters import threshold_otsu

from .common import load_config, read3d, robust_u8, write3d, write_json


# ---------------------------------------------------------------- BG1: e13 scale estimate (exact port)
def _e13_robust_u8(volume: np.ndarray) -> np.ndarray:
    values = volume.astype(np.float32, copy=False)
    p1, p99 = np.percentile(values, [1, 99])
    scaled = (values - p1) / (p99 - p1 + 1e-8)
    return (255 * np.clip(scaled, 0, 1)).astype(np.uint8)


def estimate_scale(gray: np.ndarray) -> dict[str, Any]:
    gray_u8 = _e13_robust_u8(gray)
    threshold = float(threshold_otsu(gray_u8))
    foreground = gray_u8 >= threshold
    foreground = ndi.binary_closing(foreground, structure=morphology.ball(1))
    foreground = ndi.binary_fill_holes(foreground)
    distance = ndi.distance_transform_edt(foreground)
    values = distance[foreground]
    values = values[values > 0]
    if values.size == 0:
        raise RuntimeError("Empty foreground; cannot estimate scale.")
    maximum = ndi.maximum_filter(distance, size=5)
    peak_mask = (distance == maximum) & foreground
    minimum_peak = max(2.0, float(np.percentile(values, 75)))
    peak_mask &= distance >= minimum_peak
    peak_values = distance[peak_mask]
    if peak_values.size < 10:
        radius = float(np.percentile(values, 95))
        method = "fallback_distance_p95"
    else:
        peak_values = np.sort(peak_values)[::-1][:1500]
        radius = float(np.median(peak_values))
        method = "median_local_distance_peaks"
    radius = max(radius, 3.0)
    volume = max((4.0 / 3.0) * math.pi * radius ** 3, 500.0)
    return {
        "estimated_grain_radius": radius,
        "estimated_grain_volume": volume,
        "radius_estimation_method": method,
        "n_distance_peaks_used": int(peak_values.size),
        "otsu_threshold_u8": threshold,
        "openrec_radius_e13": int(max(8, round(1.5 * radius))),
        "min_predicted_size_e13": int(max(500, round(0.10 * volume))),
        "gt_used": False,
    }


# ---------------------------------------------------------------- BG2-BG5, BG7: opening by reconstruction
def openrec_radius(cfg: dict[str, Any], scale: dict[str, Any]) -> int:
    if cfg["radius"] == "auto":
        return int(max(cfg["radius_min"], round(cfg["radius_factor"] * scale["estimated_grain_radius"])))
    return int(cfg["radius"])


def _footprint(radius: int, kind: str):
    if kind == "ball_sequence":
        return morphology.ball(radius, decomposition="sequence")
    if kind == "ball":
        return morphology.ball(radius)
    raise ValueError(f"Unknown openrec footprint {kind!r}")


def openrec_flatten(gray: np.ndarray, radius: int, cfg: dict[str, Any]) -> tuple[np.ndarray, dict[str, Any]]:
    if cfg.get("input", "raw") != "raw":
        raise ValueError("Only openrec.input=raw is implemented (BG7)")
    image = gray.astype(np.float32)
    footprint = _footprint(radius, cfg["footprint"])
    if cfg["operation"] == "opening":
        marker = morphology.erosion(image, footprint)
        background = morphology.reconstruction(marker, image, method="dilation")
    elif cfg["operation"] == "closing":
        marker = morphology.dilation(image, footprint)
        background = morphology.reconstruction(marker, image, method="erosion")
    else:
        raise ValueError(f"Unknown openrec operation {cfg['operation']!r}")

    if cfg["mode"] == "subtract_background":
        flat = image - background + float(np.median(background))
        flat = np.clip(flat, 0, None)
    elif cfg["mode"] == "reconstructed_image":
        flat = background
    else:
        raise ValueError(f"Unknown openrec mode {cfg['mode']!r}")
    p_low, p_high = cfg.get("robust_percentiles", [0.5, 99.5])
    flat_u8 = robust_u8(flat, p_low, p_high)
    stats = {
        "background_p05_p50_p95": [float(v) for v in np.percentile(background, [5, 50, 95])],
        "background_median": float(np.median(background)),
        "flat_u8_mean": float(flat_u8.mean()),
        "flat_u8_std": float(flat_u8.std()),
    }
    return flat_u8, stats


# ---------------------------------------------------------------- BG6: lambda flattening (exact port)
def flatten_volume_u8(vol_u8: np.ndarray, lambda_xy: float = 15.0, sigma_z: float = 0.0) -> np.ndarray:
    vol = vol_u8.astype(np.float32)
    if sigma_z > 0:
        bg = ndi.gaussian_filter(vol, sigma=(sigma_z, lambda_xy, lambda_xy))
    else:
        bg = np.empty_like(vol)
        for z in range(vol.shape[0]):
            bg[z] = ndi.gaussian_filter(vol[z], sigma=lambda_xy)
    bg_ref = np.median(bg[bg > 0]) if np.any(bg > 0) else np.median(bg)
    flat = np.clip(vol - bg + bg_ref, 0, None)
    return robust_u8(flat, p_low=0.5, p_high=99.5)


def lambda_tag(value: float) -> str:
    return str(float(value)).replace(".", "p")


def run(gray: np.ndarray, config: dict[str, Any], out_dir: Path, prefix: str) -> dict[str, Any]:
    out_dir = Path(out_dir)
    scale = estimate_scale(gray)
    radius = openrec_radius(config["openrec"], scale)
    openrec_u8, openrec_stats = openrec_flatten(gray, radius, config["openrec"])
    # batos_preprocess.main: robust_u8(input) (copy for uint8) then flatten with lambda
    selector = flatten_volume_u8(robust_u8(openrec_u8), config["selector"]["lambda_xy"],
                                 config["selector"]["sigma_z"])
    openrec_path = write3d(out_dir / f"{prefix}_flatten_openrec_r{radius}_8u.tif", openrec_u8)
    selector_path = write3d(
        out_dir / f"{prefix}_flatten_lambda{lambda_tag(config['selector']['lambda_xy'])}_8u.tif", selector)
    info = {
        "scale": scale,
        "openrec_radius": radius,
        "openrec": openrec_stats,
        "selector_lambda_xy": config["selector"]["lambda_xy"],
        "openrec_image": str(openrec_path),
        "selector_image": str(selector_path),
    }
    write_json(out_dir / f"{prefix}_preprocess.json", info)
    return info


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="batos-paper-v2 stage 1: preprocessing")
    parser.add_argument("gray", type=Path, help="raw gray crop (uint16)")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--config", type=Path, default=Path("configs/batos_paper_v2.yaml"))
    args = parser.parse_args(argv)
    info = run(read3d(args.gray, "gray"), load_config(args.config), args.out_dir, args.prefix)
    print(f"[preprocess] R_hat={info['scale']['estimated_grain_radius']:.3f} "
          f"openrec_r={info['openrec_radius']} -> {info['openrec_image']}, {info['selector_image']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
