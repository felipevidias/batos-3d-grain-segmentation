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
from skimage.feature import peak_local_max
from skimage.segmentation import watershed

ROOT = Path(".")
SCRIPT_DIR = ROOT / "code/segmentation/batos"
TILE_DIR = ROOT / "results_batos/full_image_pilot_tiles/EFRGP01_01_tile_center_256_halo_64"
BASE = TILE_DIR / "batos_ls_auto_prediction"

ORIGINAL_8U = BASE / "01_preprocess/pilot_tile_original_8u.tif"
FLAT = BASE / "01_preprocess/pilot_tile_flatten_lambda15p0_8u.tif"
INTERNAL_DARK = BASE / "03_background/pilot_tile_bgreject_internal_dark_candidates.tif"
GT_CENTER = TILE_DIR / "label_center_256.tif"

HALO = 64
CENTER_SIZE = 256
MIN_SIZE = 500

PYTHON = sys.executable


def run(cmd):
    print()
    print("RUN:", " ".join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), check=True)


def require(p):
    p = Path(p)
    if not p.exists():
        raise FileNotFoundError(p)
    return p


def read_labels(p):
    arr = tiff.imread(str(p))
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D, got {arr.shape}: {p}")
    return arr.astype(np.int32, copy=False)


def crop_center(arr):
    s = HALO
    e = HALO + CENTER_SIZE
    return arr[s:e, s:e, s:e]


def label_sizes(labels):
    ids, counts = np.unique(labels, return_counts=True)
    return np.array([c for i, c in zip(ids, counts) if int(i) > 0], dtype=np.float64)


def robust_auto_params(labels):
    sizes = label_sizes(labels)
    if len(sizes) == 0:
        raise RuntimeError("No labels found.")

    q10, q90 = np.percentile(sizes, [10, 90])
    clean = sizes[(sizes >= q10) & (sizes <= q90)]
    if len(clean) < 5:
        clean = sizes

    v_med = float(np.median(clean))
    r_eq = float((3.0 * v_med / (4.0 * math.pi)) ** (1.0 / 3.0))

    return {
        "grain_volume_median": v_med,
        "equivalent_radius": r_eq,
        "MERGE_MIN_VOLUME": int(round(1.5 * v_med)),
        "MIN_SPLIT_PART_SIZE": int(round(0.18 * v_med)),
        "PEAK_MIN_DISTANCE": int(max(5, round(0.35 * r_eq))),
        "PEAK_THRESHOLD": float(max(2.0, 0.18 * r_eq)),
        "MAX_PEAKS_GLOBAL_CAP": 12,
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


def split_large_labels_auto_tile(labels):
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
        "labels_topk_peaks_used": 0,
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
        max_peaks_for_label = int(min(max_peaks_cap, max(2, expected_parts + 2)))

        # mudança principal: se há picos demais, não pula; fica com os top-k pelo valor da distância
        if len(peaks) > max_peaks_for_label:
            vals = np.array([dist[z, y, x] for z, y, x in peaks], dtype=np.float32)
            order = np.argsort(vals)[::-1][:max_peaks_for_label]
            peaks = peaks[order]
            stats["labels_topk_peaks_used"] += 1

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


def colorize_labels(labels):
    labels = labels.astype(np.int64)
    rgb = np.zeros(labels.shape + (3,), dtype=np.uint8)
    ids = np.unique(labels)
    ids = ids[ids > 0]
    for lab in ids:
        r = max(int((lab * 37) % 255), 60)
        g = max(int((lab * 91) % 255), 60)
        b = max(int((lab * 53) % 255), 60)
        rgb[labels == lab] = (r, g, b)
    return rgb


def center_crop_160(arr):
    z, y, x = arr.shape
    s = 80
    return arr[z//2-s:z//2+s, y//2-s:y//2+s, x//2-s:x//2+s]


def preview(gt, batos, pred, out):
    out.mkdir(parents=True, exist_ok=True)
    gt160 = center_crop_160(gt)
    b160 = center_crop_160(batos)
    p160 = center_crop_160(pred)

    gt_rgb = colorize_labels(gt160)
    b_rgb = colorize_labels(b160)
    p_rgb = colorize_labels(p160)

    tiff.imwrite(str(out / "gt_center_160_rgb_stack.tif"), gt_rgb, photometric="rgb")
    tiff.imwrite(str(out / "batos_center_160_rgb_stack.tif"), b_rgb, photometric="rgb")
    tiff.imwrite(str(out / "batos_ls_auto_tile_center_160_rgb_stack.tif"), p_rgb, photometric="rgb")


def run_metrics(gt, pred, out):
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for th in ["0.25", "0.50", "0.75"]:
        tag = th.replace(".", "p")
        name = f"pilot_tile_BA_TOS_LS_AUTO_TILE_iou_{tag}"
        th_dir = out / f"iou_{tag}"
        cmd = [
            PYTHON, SCRIPT_DIR / "batos_instance_metrics.py",
            "--gt", str(gt),
            "--pred", str(pred),
            "--out-dir", str(th_dir),
            "--min-gt-area", "1000",
            "--ignore-gt-border",
            "--thresholds", th,
            "--name", name,
        ]
        run(cmd)
        csv = th_dir / f"{name}_instance_metrics.csv"
        row = pd.read_csv(csv).iloc[0].to_dict()
        row["threshold"] = float(th)
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(out / "multi_iou_metrics.csv", index=False)

    r025 = df[df["threshold"] == 0.25].iloc[0]
    r050 = df[df["threshold"] == 0.50].iloc[0]
    r075 = df[df["threshold"] == 0.75].iloc[0]

    summary = pd.DataFrame([{
        "method": "BA-TOS-LS-AUTO-TILE",
        "P@0.50": r050["precision"],
        "R@0.50": r050["recall"],
        "F1@0.25": r025["f1"],
        "F1@0.50": r050["f1"],
        "F1@0.75": r075["f1"],
        "overseg": r050["overseg_ratio_pred_eval_gt"],
        "tp@0.50": int(r050["tp"]),
        "fp@0.50": int(r050["fp"]),
        "fn@0.50": int(r050["fn"]),
    }])
    summary.to_csv(out / "summary.csv", index=False)
    print(summary.to_string(index=False))


def main():
    require(ORIGINAL_8U)
    require(FLAT)
    require(INTERNAL_DARK)
    require(GT_CENTER)

    for max_markers in [1400, 1800]:
        tag = f"maxmarkers_{max_markers}"
        out = TILE_DIR / f"batos_ls_auto_tile_rerun_{tag}"
        out.mkdir(parents=True, exist_ok=True)

        mk_dir = out / "04_markers"
        ws_dir = out / "05_watershed"
        ls_dir = out / "06_local_split_auto_tile"
        center_dir = out / "07_center_crop"
        metrics_dir = out / "08_metrics"
        preview_dir = out / "09_gmic_preview"

        run([
            PYTHON, SCRIPT_DIR / "batos_markers.py",
            ORIGINAL_8U,
            "--flattened", FLAT,
            "--internal-dark", INTERNAL_DARK,
            "--out-dir", mk_dir,
            "--prefix", "pilot_tile_batos",
            "--max-markers", str(max_markers),
        ])

        markers = mk_dir / "pilot_tile_batos_batos_marker_seed16_for_macke.tif"
        require(markers)

        run([
            PYTHON, SCRIPT_DIR / "batos_watershed_distance.py",
            "--image", ORIGINAL_8U,
            "--markers", markers,
            "--out-dir", ws_dir,
            "--prefix", "batos_distance",
            "--threshold-offset", "-10",
            "--close-radius", "2",
            "--force-markers-foreground",
            "--force-radius", "1",
        ])

        raw_labels = ws_dir / "batos_distance_labels.tif"
        require(raw_labels)

        run([
            PYTHON, SCRIPT_DIR / "filter_predicted_labels_by_size_labels_only.py",
            "--labels", raw_labels,
            "--out-dir", ws_dir,
            "--prefix", "batos_distance_filtered",
            "--min-size", str(MIN_SIZE),
            "--relabel",
        ])

        batos_labels = ws_dir / f"batos_distance_filtered_min{MIN_SIZE}_labels.tif"
        require(batos_labels)

        labels = read_labels(batos_labels)
        split, stats, params = split_large_labels_auto_tile(labels)

        ls_dir.mkdir(parents=True, exist_ok=True)
        pred_proc = ls_dir / "batos_ls_auto_tile_pred_proc_384.tif"
        tiff.imwrite(str(pred_proc), split.astype(np.uint16))

        with open(ls_dir / "params.json", "w") as f:
            json.dump(params, f, indent=2)
        with open(ls_dir / "stats.json", "w") as f:
            json.dump(stats, f, indent=2)

        batos_center = crop_center(labels)
        pred_center = crop_center(split)
        gt_center = read_labels(GT_CENTER)

        center_dir.mkdir(parents=True, exist_ok=True)
        batos_center_path = center_dir / "batos_pred_center_256.tif"
        pred_center_path = center_dir / "batos_ls_auto_tile_pred_center_256.tif"
        gt_center_path = center_dir / "gt_label_center_256.tif"

        tiff.imwrite(str(batos_center_path), batos_center.astype(np.uint16))
        tiff.imwrite(str(pred_center_path), pred_center.astype(np.uint16))
        tiff.imwrite(str(gt_center_path), gt_center.astype(np.uint16))

        print()
        print("===", tag, "params ===")
        print(json.dumps(params, indent=2))
        print("===", tag, "stats ===")
        print(json.dumps(stats, indent=2))

        run_metrics(gt_center_path, pred_center_path, metrics_dir)
        preview(gt_center, batos_center, pred_center, preview_dir)

        print()
        print("Open prediction alone:")
        print(f"gmic {preview_dir / 'batos_ls_auto_tile_center_160_rgb_stack.tif'} -append z -display")
        print("Open GT alone:")
        print(f"gmic {preview_dir / 'gt_center_160_rgb_stack.tif'} -append z -display")


if __name__ == "__main__":
    main()
