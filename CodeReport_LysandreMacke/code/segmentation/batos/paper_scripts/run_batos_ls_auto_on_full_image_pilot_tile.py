#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys
import json
import math

import numpy as np
import pandas as pd
import tifffile as tiff
from scipy import ndimage as ndi

try:
    from skimage.feature import peak_local_max
    from skimage.segmentation import watershed
except Exception as e:
    raise SystemExit(
        "scikit-image is required. Install with: python -m pip install scikit-image\n"
        f"Original error: {e}"
    )


ROOT = Path(".")
SCRIPT_DIR = ROOT / "code/segmentation/batos"

TILE_DIR = ROOT / "results_batos/full_image_pilot_tiles/EFRGP01_01_tile_center_256_halo_64"

GRAY_PROC = TILE_DIR / "gray_proc_384_with_halo.tif"
GT_PROC = TILE_DIR / "label_proc_384_with_halo.tif"
GRAY_CENTER = TILE_DIR / "gray_center_256.tif"
GT_CENTER = TILE_DIR / "label_center_256.tif"

OUT = TILE_DIR / "batos_ls_auto_prediction"
OUT.mkdir(parents=True, exist_ok=True)

HALO = 64
CENTER_SIZE = 256
MIN_SIZE = 500

PYTHON = sys.executable


def run(cmd):
    print()
    print("RUN:", " ".join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), check=True)


def find_one(pattern):
    hits = sorted(Path(".").glob(str(pattern)))
    if not hits:
        return None
    return hits[0]


def require(path, label):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def read_labels(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D label volume, got {arr.shape}: {path}")
    return arr.astype(np.int32, copy=False)


def label_sizes(labels):
    ids, counts = np.unique(labels, return_counts=True)
    return np.array([c for i, c in zip(ids, counts) if int(i) > 0], dtype=np.float64)


def robust_auto_params(labels):
    sizes = label_sizes(labels)

    if len(sizes) == 0:
        raise RuntimeError("No labels found to estimate automatic parameters.")

    q10, q90 = np.percentile(sizes, [10, 90])
    clean = sizes[(sizes >= q10) & (sizes <= q90)]

    if len(clean) < 5:
        clean = sizes

    v_med = float(np.median(clean))
    v_q25 = float(np.percentile(clean, 25))
    v_q75 = float(np.percentile(clean, 75))
    r_eq = float((3.0 * v_med / (4.0 * math.pi)) ** (1.0 / 3.0))

    return {
        "grain_volume_median": v_med,
        "grain_volume_q25": v_q25,
        "grain_volume_q75": v_q75,
        "equivalent_radius": r_eq,
        "MERGE_MIN_VOLUME": int(round(1.8 * v_med)),
        "MIN_SPLIT_PART_SIZE": int(round(0.25 * v_med)),
        "PEAK_MIN_DISTANCE": int(max(3, round(0.45 * r_eq))),
        "PEAK_THRESHOLD": float(max(2.0, 0.25 * r_eq)),
        "MAX_PEAKS_GLOBAL_CAP": 6,
    }


def bbox_from_mask(mask, margin=2):
    coords = np.argwhere(mask)
    if coords.size == 0:
        return None

    z0, y0, x0 = coords.min(axis=0)
    z1, y1, x1 = coords.max(axis=0) + 1

    z0 = max(0, z0 - margin)
    y0 = max(0, y0 - margin)
    x0 = max(0, x0 - margin)

    z1 = min(mask.shape[0], z1 + margin)
    y1 = min(mask.shape[1], y1 + margin)
    x1 = min(mask.shape[2], x1 + margin)

    return slice(z0, z1), slice(y0, y1), slice(x0, x1)


def split_large_labels_auto(labels):
    params = robust_auto_params(labels)

    merge_min_volume = params["MERGE_MIN_VOLUME"]
    min_split_part_size = params["MIN_SPLIT_PART_SIZE"]
    peak_min_distance = params["PEAK_MIN_DISTANCE"]
    peak_threshold = params["PEAK_THRESHOLD"]
    max_peaks_cap = params["MAX_PEAKS_GLOBAL_CAP"]
    v_med = params["grain_volume_median"]

    out = np.zeros_like(labels, dtype=np.int32)
    next_id = 1

    stats = {
        "original_labels": 0,
        "new_labels": 0,
        "candidate_large_labels": 0,
        "labels_split": 0,
        "labels_kept_unsplit": 0,
        "labels_skipped_too_many_peaks": 0,
        "small_split_parts_removed": 0,
    }

    ids, counts = np.unique(labels, return_counts=True)
    lab_sizes = {int(i): int(c) for i, c in zip(ids, counts) if int(i) > 0}

    for lab in sorted(lab_sizes):
        stats["original_labels"] += 1
        size = lab_sizes[lab]
        mask = labels == lab

        if size < merge_min_volume:
            out[mask] = next_id
            next_id += 1
            stats["new_labels"] += 1
            stats["labels_kept_unsplit"] += 1
            continue

        stats["candidate_large_labels"] += 1

        bb = bbox_from_mask(mask, margin=2)
        if bb is None:
            continue

        local_mask = mask[bb]
        dist = ndi.distance_transform_edt(local_mask)

        peaks = peak_local_max(
            dist,
            labels=local_mask.astype(np.uint8),
            min_distance=peak_min_distance,
            threshold_abs=peak_threshold,
            exclude_border=False,
        )

        if len(peaks) < 2:
            out[mask] = next_id
            next_id += 1
            stats["new_labels"] += 1
            stats["labels_kept_unsplit"] += 1
            continue

        expected_parts = int(math.ceil(size / max(v_med, 1.0)))
        max_peaks_for_label = int(min(max_peaks_cap, max(2, expected_parts + 1)))

        if len(peaks) > max_peaks_for_label:
            out[mask] = next_id
            next_id += 1
            stats["new_labels"] += 1
            stats["labels_skipped_too_many_peaks"] += 1
            continue

        markers = np.zeros(local_mask.shape, dtype=np.int32)

        for k, (z, y, x) in enumerate(peaks, start=1):
            markers[z, y, x] = k

        split = watershed(-dist, markers, mask=local_mask)

        kept_parts = []

        for sid in range(1, int(split.max()) + 1):
            part = split == sid
            part_size = int(part.sum())

            if part_size >= min_split_part_size:
                kept_parts.append(part)
            else:
                stats["small_split_parts_removed"] += 1

        if len(kept_parts) < 2:
            out[mask] = next_id
            next_id += 1
            stats["new_labels"] += 1
            stats["labels_kept_unsplit"] += 1
            continue

        local_out = np.zeros(local_mask.shape, dtype=np.int32)

        for part in kept_parts:
            local_out[part] = next_id
            next_id += 1
            stats["new_labels"] += 1

        out_bb = out[bb]
        out_bb[local_out > 0] = local_out[local_out > 0]
        out[bb] = out_bb

        stats["labels_split"] += 1

    return out, stats, params


def crop_center_from_proc(arr):
    start = HALO
    end = HALO + CENTER_SIZE
    return arr[start:end, start:end, start:end]


def colorize_labels(labels):
    labels = labels.astype(np.int64)
    rgb = np.zeros(labels.shape + (3,), dtype=np.uint8)

    ids = np.unique(labels)
    ids = ids[ids > 0]

    for lab in ids:
        r = int((lab * 37) % 255)
        g = int((lab * 91) % 255)
        b = int((lab * 53) % 255)

        r = max(r, 60)
        g = max(g, 60)
        b = max(b, 60)

        rgb[labels == lab] = (r, g, b)

    return rgb


def center_crop(arr, size=160):
    z, y, x = arr.shape
    z0 = max(0, z // 2 - size // 2)
    y0 = max(0, y // 2 - size // 2)
    x0 = max(0, x // 2 - size // 2)
    return arr[z0:z0+size, y0:y0+size, x0:x0+size]


def make_gmic_color_stack(gt_center, pred_center, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)

    gt160 = center_crop(gt_center, 160)
    pred160 = center_crop(pred_center, 160)

    gt_rgb = colorize_labels(gt160)
    pred_rgb = colorize_labels(pred160)

    tiff.imwrite(str(out_dir / "gt_center_160_rgb_stack.tif"), gt_rgb, photometric="rgb")
    tiff.imwrite(str(out_dir / "pred_batos_ls_auto_center_160_rgb_stack.tif"), pred_rgb, photometric="rgb")

    gap = np.ones((gt_rgb.shape[0], gt_rgb.shape[1], 8, 3), dtype=np.uint8) * 255
    side = np.concatenate([gt_rgb, gap, pred_rgb], axis=2)

    tiff.imwrite(str(out_dir / "gt_vs_batos_ls_auto_center_160_rgb_stack.tif"), side, photometric="rgb")

    print()
    print("G'MIC preview:")
    print(out_dir / "gt_center_160_rgb_stack.tif")
    print(out_dir / "pred_batos_ls_auto_center_160_rgb_stack.tif")
    print(out_dir / "gt_vs_batos_ls_auto_center_160_rgb_stack.tif")


def main():
    require(GRAY_PROC, "GRAY_PROC")
    require(GT_PROC, "GT_PROC")
    require(GRAY_CENTER, "GRAY_CENTER")
    require(GT_CENTER, "GT_CENTER")

    print("=== BA-TOS-LS-AUTO full-image pilot tile ===")
    print("GRAY_PROC :", GRAY_PROC)
    print("GT_PROC   :", GT_PROC)
    print("GRAY_CENTER:", GRAY_CENTER)
    print("GT_CENTER :", GT_CENTER)
    print("OUT       :", OUT)

    # 01 preprocess
    pp_dir = OUT / "01_preprocess"
    run([
        PYTHON, SCRIPT_DIR / "batos_preprocess.py",
        GRAY_PROC,
        "--out-dir", pp_dir,
        "--prefix", "pilot_tile",
        "--lambda-xy", "15",
        "--sigma-z", "0",
    ])

    original_8u = pp_dir / "pilot_tile_original_8u.tif"
    flat_candidates = sorted(pp_dir.glob("pilot_tile_flatten_lambda*_8u.tif"))
    if not flat_candidates:
        raise FileNotFoundError(f"No flattened file found in {pp_dir}")
    flat = flat_candidates[0]
    print("selected flat:", flat)

    require(original_8u, "original_8u")
    require(flat, "flat")

    # 02 ToS seed
    tos_dir = OUT / "02_tos_seed"
    run([
        PYTHON, SCRIPT_DIR / "batos_tos_seed.py",
        flat,
        "--out-dir", tos_dir,
    ])

    seed_candidates = sorted(tos_dir.glob("*.tif"))
    if not seed_candidates:
        raise RuntimeError(f"No ToS seed tif generated in {tos_dir}")

    # tenta achar seed binária/core mais provável
    preferred = [
        p for p in seed_candidates
        if "seed" in p.name.lower() or "core" in p.name.lower() or "mask" in p.name.lower()
    ]
    raw_seed = preferred[0] if preferred else seed_candidates[0]

    print("selected raw ToS seed:", raw_seed)

    # 03 boundary-aware background rejection
    bg_dir = OUT / "03_background"
    run([
        PYTHON, SCRIPT_DIR / "batos_background.py",
        original_8u,
        "--flattened", flat,
        "--seed", raw_seed,
        "--out-dir", bg_dir,
        "--prefix", "pilot_tile_bgreject",
    ])

    internal_candidates = sorted(bg_dir.glob("*internal*dark*candidates*.tif"))
    if not internal_candidates:
        internal_candidates = sorted(bg_dir.glob("*internal*.tif"))

    if not internal_candidates:
        raise RuntimeError(f"No internal dark candidate tif found in {bg_dir}")

    internal_dark = internal_candidates[0]
    print("selected internal dark:", internal_dark)

    # 04 markers
    mk_dir = OUT / "04_markers"
    run([
        PYTHON, SCRIPT_DIR / "batos_markers.py",
        original_8u,
        "--flattened", flat,
        "--internal-dark", internal_dark,
        "--out-dir", mk_dir,
        "--prefix", "pilot_tile_batos",
        "--max-markers", "900",
    ])

    marker_candidates = sorted(mk_dir.glob("*marker*seed16*.tif"))
    if not marker_candidates:
        marker_candidates = sorted(mk_dir.glob("*centroid*.tif"))
    if not marker_candidates:
        marker_candidates = sorted(mk_dir.glob("*.tif"))

    if not marker_candidates:
        raise RuntimeError(f"No marker tif found in {mk_dir}")

    markers = marker_candidates[0]
    print("selected markers:", markers)

    # 05 watershed
    ws_dir = OUT / "05_watershed"
    run([
        PYTHON, SCRIPT_DIR / "batos_watershed_distance.py",
        "--image", original_8u,
        "--markers", markers,
        "--out-dir", ws_dir,
        "--prefix", "batos_distance",
        "--threshold-offset", "-10",
        "--close-radius", "2",
        "--force-markers-foreground",
        "--force-radius", "1",
    ])

    raw_labels = ws_dir / "batos_distance_labels.tif"
    require(raw_labels, "raw watershed labels")

    # 06 size filter
    run([
        PYTHON, SCRIPT_DIR / "filter_predicted_labels_by_size_labels_only.py",
        "--labels", raw_labels,
        "--out-dir", ws_dir,
        "--prefix", "batos_distance_filtered",
        "--min-size", str(MIN_SIZE),
        "--relabel",
    ])

    filtered_labels = ws_dir / f"batos_distance_filtered_min{MIN_SIZE}_labels.tif"
    require(filtered_labels, "filtered BA-TOS labels")

    # 07 local split auto
    ls_dir = OUT / "06_local_split_auto"
    ls_dir.mkdir(parents=True, exist_ok=True)

    batos_labels = read_labels(filtered_labels)
    split_labels, stats, params = split_large_labels_auto(batos_labels)

    auto_proc_path = ls_dir / "batos_ls_auto_pred_proc_384.tif"
    tiff.imwrite(str(auto_proc_path), split_labels.astype(np.uint16))

    with open(ls_dir / "batos_ls_auto_params.json", "w") as f:
        json.dump(params, f, indent=2)

    with open(ls_dir / "batos_ls_auto_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    print()
    print("BA-TOS-LS-AUTO params:")
    print(json.dumps(params, indent=2))
    print("BA-TOS-LS-AUTO stats:")
    print(json.dumps(stats, indent=2))
    print("wrote:", auto_proc_path)

    # 08 crop center
    center_dir = OUT / "07_center_crop"
    center_dir.mkdir(parents=True, exist_ok=True)

    pred_center = crop_center_from_proc(split_labels)
    batos_center = crop_center_from_proc(batos_labels)
    gt_center = read_labels(GT_CENTER)

    pred_center_path = center_dir / "batos_ls_auto_pred_center_256.tif"
    batos_center_path = center_dir / "batos_pred_center_256.tif"
    gt_center_copy = center_dir / "gt_label_center_256.tif"

    tiff.imwrite(str(pred_center_path), pred_center.astype(np.uint16))
    tiff.imwrite(str(batos_center_path), batos_center.astype(np.uint16))
    tiff.imwrite(str(gt_center_copy), gt_center.astype(np.uint16))

    print("wrote:", batos_center_path)
    print("wrote:", pred_center_path)
    print("wrote:", gt_center_copy)

    # 09 metrics center
    metrics_dir = OUT / "08_metrics"
    run([
        PYTHON, SCRIPT_DIR / "batos_instance_metrics.py",
        "--gt", gt_center_copy,
        "--pred", pred_center_path,
        "--out-dir", metrics_dir,
        "--min-gt-area", "1000",
        "--ignore-gt-border",
        "--thresholds", "0.25", "0.50", "0.75",
        "--name", "pilot_tile_center_256_BA_TOS_LS_AUTO",
    ])

    # 10 gmic preview
    make_gmic_color_stack(
        gt_center=gt_center,
        pred_center=pred_center,
        out_dir=OUT / "09_gmic_preview",
    )

    # manifest
    manifest = {
        "gray_proc": str(GRAY_PROC),
        "gt_proc": str(GT_PROC),
        "gray_center": str(GRAY_CENTER),
        "gt_center": str(GT_CENTER),
        "original_8u": str(original_8u),
        "flat": str(flat),
        "raw_seed": str(raw_seed),
        "internal_dark": str(internal_dark),
        "markers": str(markers),
        "raw_labels": str(raw_labels),
        "filtered_batos_labels": str(filtered_labels),
        "batos_ls_auto_proc": str(auto_proc_path),
        "batos_center": str(batos_center_path),
        "batos_ls_auto_center": str(pred_center_path),
        "gt_center_copy": str(gt_center_copy),
        "metrics_dir": str(metrics_dir),
        "gmic_preview": str(OUT / "09_gmic_preview"),
        "halo": HALO,
        "center_size": CENTER_SIZE,
        "min_size": MIN_SIZE,
        "auto_params": params,
        "auto_stats": stats,
    }

    with open(OUT / "pilot_tile_batos_ls_auto_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    print()
    print("DONE.")
    print("Manifest:", OUT / "pilot_tile_batos_ls_auto_manifest.json")


if __name__ == "__main__":
    main()
