#!/usr/bin/env python3
from pathlib import Path
import csv
import json
import math
import subprocess
import sys

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
MANIFEST = ROOT / "results_batos/paper_tables/final_revision_clean/manual_clean_manifest_multi_iou.csv"
METRICS = ROOT / "code/segmentation/batos/batos_instance_metrics.py"

OUT = ROOT / "results_batos/paper_tables/batos_local_split_auto"
OUT.mkdir(parents=True, exist_ok=True)

METHOD = "BA-TOS"
THRESHOLDS = ["0.25", "0.50", "0.75"]


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

    params = {
        "grain_volume_median": v_med,
        "grain_volume_q25": v_q25,
        "grain_volume_q75": v_q75,
        "equivalent_radius": r_eq,

        # escala automática
        "MERGE_MIN_VOLUME": int(round(1.8 * v_med)),
        "MIN_SPLIT_PART_SIZE": int(round(0.25 * v_med)),
        "PEAK_MIN_DISTANCE": int(max(3, round(0.45 * r_eq))),
        "PEAK_THRESHOLD": float(max(2.0, 0.25 * r_eq)),
        "MAX_PEAKS_GLOBAL_CAP": 6,
    }

    return params


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

        # número máximo esperado de partes baseado no volume do label
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


def load_manifest():
    rows = list(csv.DictReader(open(MANIFEST)))
    clean = []

    for r in rows:
        method = r.get("method") or r.get("variant")

        if method != METHOD:
            continue

        if r.get("gt_found", "1") == "0" or r.get("pred_found", "1") == "0":
            continue

        clean.append(r)

    return clean


def run_metrics(gt, pred, out_dir, name, threshold):
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, str(METRICS),
        "--gt", str(gt),
        "--pred", str(pred),
        "--out-dir", str(out_dir),
        "--min-gt-area", "1000",
        "--ignore-gt-border",
        "--thresholds", threshold,
        "--name", name,
    ]

    print("RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    metric_file = out_dir / f"{name}_instance_metrics.csv"

    if not metric_file.exists():
        raise FileNotFoundError(f"Metric file not found: {metric_file}")

    return pd.read_csv(metric_file).iloc[0].to_dict()


def main():
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Manifest not found: {MANIFEST}")

    if not METRICS.exists():
        raise FileNotFoundError(f"Metrics script not found: {METRICS}")

    rows = load_manifest()

    print("BA-TOS rows:", len(rows))
    print("Output:", OUT)

    manifest_rows = []
    metrics_rows = []

    for r in rows:
        scan = r["scan"]
        crop_id = r["crop_id"]
        gt = Path(r["gt"])
        pred = Path(r["pred"])

        labels = read_labels(pred)
        split_labels, stats, params = split_large_labels_auto(labels)

        crop_out = OUT / scan / crop_id
        crop_out.mkdir(parents=True, exist_ok=True)

        split_path = crop_out / f"{scan}_{crop_id}_batos_local_split_AUTO_labels.tif"
        tiff.imwrite(str(split_path), split_labels.astype(np.int32))

        param_path = crop_out / f"{scan}_{crop_id}_auto_params.json"
        with open(param_path, "w") as f:
            json.dump(params, f, indent=2)

        manifest_row = {
            "scan": scan,
            "crop_id": crop_id,
            "gt": str(gt),
            "original_pred": str(pred),
            "split_pred": str(split_path),
            "auto_params": str(param_path),
            **params,
            **stats,
        }

        manifest_rows.append(manifest_row)

        for th in THRESHOLDS:
            tag = th.replace(".", "p")
            name = f"{scan}_{crop_id}_BA_TOS_LS_AUTO_iou_{tag}"
            metric_out = crop_out / f"metrics_iou_{tag}"

            m = run_metrics(
                gt=gt,
                pred=split_path,
                out_dir=metric_out,
                name=name,
                threshold=th,
            )

            metrics_rows.append({
                "method": "BA-TOS-LS-AUTO",
                "scan": scan,
                "crop_id": crop_id,
                "threshold": float(th),
                **m,
            })

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_csv = OUT / "local_split_auto_manifest.csv"
    manifest_df.to_csv(manifest_csv, index=False)

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_csv = OUT / "local_split_auto_crop_metrics.csv"
    metrics_df.to_csv(metrics_csv, index=False)

    summary_rows = []

    for th in [0.25, 0.50, 0.75]:
        sub = metrics_df[metrics_df["threshold"] == th].copy()

        summary_rows.append({
            "method": "BA-TOS-LS-AUTO",
            "threshold": th,
            "mean_precision": sub["precision"].astype(float).mean(),
            "mean_recall": sub["recall"].astype(float).mean(),
            "mean_f1": sub["f1"].astype(float).mean(),
            "mean_overseg": sub["overseg_ratio_pred_eval_gt"].astype(float).mean(),
            "total_fp": int(sub["fp"].astype(int).sum()),
            "total_fn": int(sub["fn"].astype(int).sum()),
            "n_crops": int(len(sub)),
        })

    by_threshold = pd.DataFrame(summary_rows)
    by_threshold_csv = OUT / "local_split_auto_by_threshold.csv"
    by_threshold.to_csv(by_threshold_csv, index=False)

    row025 = by_threshold[by_threshold["threshold"] == 0.25].iloc[0]
    row050 = by_threshold[by_threshold["threshold"] == 0.50].iloc[0]
    row075 = by_threshold[by_threshold["threshold"] == 0.75].iloc[0]

    final = pd.DataFrame([{
        "method": "BA-TOS-LS-AUTO",
        "mean_precision": row050["mean_precision"],
        "mean_recall": row050["mean_recall"],
        "F1@0.25": row025["mean_f1"],
        "F1@0.50": row050["mean_f1"],
        "F1@0.75": row075["mean_f1"],
        "mean_overseg": row050["mean_overseg"],
        "total_fp": int(row050["total_fp"]),
        "total_fn": int(row050["total_fn"]),
        "n_crops": int(row050["n_crops"]),
    }])

    final_csv = OUT / "local_split_auto_summary.csv"
    final.to_csv(final_csv, index=False)

    print()
    print("=== BA-TOS-LS-AUTO summary ===")
    print(final.to_string(index=False))
    print()
    print("wrote:", manifest_csv)
    print("wrote:", metrics_csv)
    print("wrote:", by_threshold_csv)
    print("wrote:", final_csv)


if __name__ == "__main__":
    main()
