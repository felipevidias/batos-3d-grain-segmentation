#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import tifffile as tiff


ROOT = Path(".")

MANIFEST = ROOT / "results_batos/paper_tables/batos_local_split_auto/local_split_auto_manifest.csv"

OUT = ROOT / "results_batos/paper_tables/grain_centered_analysis_batos_ls_auto"
OUT.mkdir(parents=True, exist_ok=True)

METHOD = "BA-TOS-LS-AUTO"

MIN_GT_VOLUME = 1000
IOU_050 = 0.50
IOU_075 = 0.75


def read_labels(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D label volume, got shape={arr.shape}: {path}")
    return arr.astype(np.int32, copy=False)


def labels_touching_border(labels):
    border_ids = set()

    faces = [
        labels[0, :, :],
        labels[-1, :, :],
        labels[:, 0, :],
        labels[:, -1, :],
        labels[:, :, 0],
        labels[:, :, -1],
    ]

    for face in faces:
        ids = np.unique(face)
        for lab in ids:
            lab = int(lab)
            if lab > 0:
                border_ids.add(lab)

    return border_ids


def label_stats(labels):
    ids, counts = np.unique(labels, return_counts=True)

    volumes = {}
    centroids = {}

    for lab, count in zip(ids, counts):
        lab = int(lab)
        if lab <= 0:
            continue

        mask = labels == lab
        coords = np.argwhere(mask)

        if coords.size == 0:
            continue

        volumes[lab] = int(count)
        centroids[lab] = coords.mean(axis=0)

    return volumes, centroids


def best_match_for_gt(gt_lab, gt, pred, gt_volume, pred_volumes, gt_centroid, pred_centroids):
    gt_mask = gt == gt_lab

    pred_inside = pred[gt_mask]
    candidate_ids, intersections = np.unique(pred_inside, return_counts=True)

    best_pred = 0
    best_iou = 0.0
    best_intersection = 0

    for pred_lab, inter in zip(candidate_ids, intersections):
        pred_lab = int(pred_lab)
        inter = int(inter)

        if pred_lab <= 0:
            continue

        pred_volume = pred_volumes.get(pred_lab, 0)

        if pred_volume <= 0:
            continue

        union = gt_volume + pred_volume - inter

        if union <= 0:
            continue

        iou = inter / union

        if iou > best_iou:
            best_iou = float(iou)
            best_pred = pred_lab
            best_intersection = inter

    if best_pred > 0:
        pred_volume = pred_volumes[best_pred]
        pred_centroid = pred_centroids[best_pred]

        centroid_error = float(np.linalg.norm(gt_centroid - pred_centroid))
        relative_volume_error = float(abs(pred_volume - gt_volume) / max(gt_volume, 1))
    else:
        pred_volume = 0
        centroid_error = np.nan
        relative_volume_error = np.nan

    return {
        "best_pred_label": int(best_pred),
        "best_iou": float(best_iou),
        "intersection": int(best_intersection),
        "gt_volume": int(gt_volume),
        "pred_volume": int(pred_volume),
        "recovered_iou050": int(best_iou >= IOU_050),
        "recovered_iou075": int(best_iou >= IOU_075),
        "centroid_error": centroid_error,
        "relative_volume_error": relative_volume_error,
    }


def analyze_one_crop(row):
    scan = str(row["scan"])
    crop_id = str(row["crop_id"])

    gt_path = Path(str(row["gt"]))
    pred_path = Path(str(row["split_pred"]))

    if not gt_path.exists():
        raise FileNotFoundError(f"GT not found: {gt_path}")

    if not pred_path.exists():
        raise FileNotFoundError(f"Prediction not found: {pred_path}")

    gt = read_labels(gt_path)
    pred = read_labels(pred_path)

    if gt.shape != pred.shape:
        raise RuntimeError(
            f"Shape mismatch for {scan} {crop_id}: gt={gt.shape}, pred={pred.shape}"
        )

    gt_volumes, gt_centroids = label_stats(gt)
    pred_volumes, pred_centroids = label_stats(pred)

    border_gt = labels_touching_border(gt)

    rows = []

    for gt_lab in sorted(gt_volumes):
        gt_volume = gt_volumes[gt_lab]

        if gt_volume < MIN_GT_VOLUME:
            continue

        if gt_lab in border_gt:
            continue

        match = best_match_for_gt(
            gt_lab=gt_lab,
            gt=gt,
            pred=pred,
            gt_volume=gt_volume,
            pred_volumes=pred_volumes,
            gt_centroid=gt_centroids[gt_lab],
            pred_centroids=pred_centroids,
        )

        rows.append({
            "method": METHOD,
            "scan": scan,
            "crop_id": crop_id,
            "gt_label": int(gt_lab),
            **match,
        })

    return rows


def main():
    if not MANIFEST.exists():
        raise FileNotFoundError(
            f"Manifest not found: {MANIFEST}\n"
            "Run this script from the project root."
        )

    manifest = pd.read_csv(MANIFEST)

    all_rows = []

    for _, row in manifest.iterrows():
        scan = str(row["scan"])
        crop_id = str(row["crop_id"])

        print(f"[RUN] {scan} {crop_id}")

        crop_rows = analyze_one_crop(row)
        all_rows.extend(crop_rows)

        print(f"      valid grains: {len(crop_rows)}")

    grain_df = pd.DataFrame(all_rows)

    grain_csv = OUT / "grain_centered_per_grain_batos_ls_auto.csv"
    grain_df.to_csv(grain_csv, index=False)

    if len(grain_df) == 0:
        raise RuntimeError("No valid grains found.")

    summary = pd.DataFrame([{
        "method": METHOD,
        "n_grains": int(len(grain_df)),
        "recall_iou050": float(grain_df["recovered_iou050"].mean()),
        "recall_iou075": float(grain_df["recovered_iou075"].mean()),
        "mean_best_iou": float(grain_df["best_iou"].mean()),
        "median_best_iou": float(grain_df["best_iou"].median()),
        "mean_centroid_error": float(grain_df["centroid_error"].mean(skipna=True)),
        "median_centroid_error": float(grain_df["centroid_error"].median(skipna=True)),
        "mean_relative_volume_error": float(grain_df["relative_volume_error"].mean(skipna=True)),
        "median_relative_volume_error": float(grain_df["relative_volume_error"].median(skipna=True)),
    }])

    summary_csv = OUT / "grain_centered_summary_by_method.csv"
    summary.to_csv(summary_csv, index=False)

    crop_summary = (
        grain_df
        .groupby(["method", "scan", "crop_id"], as_index=False)
        .agg(
            n_grains=("gt_label", "count"),
            recall_iou050=("recovered_iou050", "mean"),
            recall_iou075=("recovered_iou075", "mean"),
            mean_best_iou=("best_iou", "mean"),
            median_best_iou=("best_iou", "median"),
            mean_centroid_error=("centroid_error", "mean"),
            median_centroid_error=("centroid_error", "median"),
            mean_relative_volume_error=("relative_volume_error", "mean"),
            median_relative_volume_error=("relative_volume_error", "median"),
        )
    )

    crop_summary_csv = OUT / "grain_centered_summary_by_crop.csv"
    crop_summary.to_csv(crop_summary_csv, index=False)

    print()
    print("=== Grain-centered summary: BA-TOS-LS-AUTO ===")
    print(summary.to_string(index=False))
    print()
    print("wrote:", grain_csv)
    print("wrote:", summary_csv)
    print("wrote:", crop_summary_csv)


if __name__ == "__main__":
    main()
