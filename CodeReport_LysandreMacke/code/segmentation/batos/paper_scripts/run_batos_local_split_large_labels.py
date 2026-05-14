#!/usr/bin/env python3
from pathlib import Path
import csv
import os
import subprocess
import numpy as np
import pandas as pd
import tifffile as tiff
from scipy import ndimage as ndi

try:
    from skimage.feature import peak_local_max
    from skimage.segmentation import watershed
except Exception as e:
    raise SystemExit(
        "scikit-image is required for this script. Install with: python -m pip install scikit-image\n"
        f"Original error: {e}"
    )

ROOT = Path(".")
MANIFEST = ROOT / "results_batos/paper_tables/final_revision_clean/manual_clean_manifest_multi_iou.csv"
OUT = ROOT / "results_batos/paper_tables/batos_local_split_large_labels"
OUT.mkdir(parents=True, exist_ok=True)

PY = "/home/felipe/Pesquisa-Grain_Seg/Filtro_Conexo/venv/bin/python3"
METRICS = ROOT / "code/segmentation/batos/batos_instance_metrics.py"

METHOD = "BA-TOS"

# Valores conservadores. Dá para ajustar por variável de ambiente.
MERGE_MIN_VOLUME = int(os.environ.get("BTS_MERGE_MIN_VOLUME", "20000"))
PEAK_MIN_DISTANCE = int(os.environ.get("BTS_PEAK_MIN_DISTANCE", "12"))
PEAK_THRESHOLD = float(os.environ.get("BTS_PEAK_THRESHOLD", "6.0"))
MAX_PEAKS = int(os.environ.get("BTS_MAX_PEAKS", "4"))
MIN_SPLIT_PART_SIZE = int(os.environ.get("BTS_MIN_SPLIT_PART_SIZE", "3000"))


def read_labels(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D label volume, got {arr.shape}: {path}")
    return arr.astype(np.int32, copy=False)


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


def split_large_labels(labels):
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

        # Labels pequenos/médios são mantidos.
        if size < MERGE_MIN_VOLUME:
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

        # Distance transform dentro do label predito.
        dist = ndi.distance_transform_edt(local_mask)

        peaks = peak_local_max(
            dist,
            labels=local_mask.astype(np.uint8),
            min_distance=PEAK_MIN_DISTANCE,
            threshold_abs=PEAK_THRESHOLD,
            exclude_border=False,
        )

        # Se não há múltiplos picos confiáveis, mantém.
        if len(peaks) < 2:
            out[mask] = next_id
            next_id += 1
            stats["new_labels"] += 1
            stats["labels_kept_unsplit"] += 1
            continue

        # Se há picos demais, é arriscado fragmentar; mantém nesta versão conservadora.
        if len(peaks) > MAX_PEAKS:
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
            if part_size >= MIN_SPLIT_PART_SIZE:
                kept_parts.append(part)
            else:
                stats["small_split_parts_removed"] += 1

        # Se, após filtro, não sobram 2 partes grandes, mantém o label original.
        if len(kept_parts) < 2:
            out[mask] = next_id
            next_id += 1
            stats["new_labels"] += 1
            stats["labels_kept_unsplit"] += 1
            continue

        # Aplica split no volume global.
        local_out = np.zeros(local_mask.shape, dtype=np.int32)
        for part in kept_parts:
            local_out[part] = next_id
            next_id += 1
            stats["new_labels"] += 1

        out_bb = out[bb]
        out_bb[local_out > 0] = local_out[local_out > 0]
        out[bb] = out_bb

        stats["labels_split"] += 1

    return out, stats


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


def run_metrics(gt, pred, out_dir, name):
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        PY, str(METRICS),
        "--gt", str(gt),
        "--pred", str(pred),
        "--out-dir", str(out_dir),
        "--min-gt-area", "1000",
        "--ignore-gt-border",
        "--thresholds", "0.50",
        "--name", name,
    ]
    print("RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    print("=== local split config ===")
    print("MERGE_MIN_VOLUME     =", MERGE_MIN_VOLUME)
    print("PEAK_MIN_DISTANCE    =", PEAK_MIN_DISTANCE)
    print("PEAK_THRESHOLD       =", PEAK_THRESHOLD)
    print("MAX_PEAKS            =", MAX_PEAKS)
    print("MIN_SPLIT_PART_SIZE  =", MIN_SPLIT_PART_SIZE)
    print()

    rows = load_manifest()
    print("BA-TOS rows:", len(rows))

    manifest_rows = []

    for r in rows:
        scan = r["scan"]
        crop_id = r["crop_id"]
        gt = Path(r["gt"])
        pred = Path(r["pred"])

        labels = read_labels(pred)
        split_labels, stats = split_large_labels(labels)

        crop_out = OUT / scan / crop_id
        crop_out.mkdir(parents=True, exist_ok=True)

        split_path = crop_out / f"{scan}_{crop_id}_batos_local_split_labels.tif"
        tiff.imwrite(str(split_path), split_labels.astype(np.int32))

        row = {
            "scan": scan,
            "crop_id": crop_id,
            "gt": str(gt),
            "original_pred": str(pred),
            "split_pred": str(split_path),
            **stats,
        }
        manifest_rows.append(row)

        run_metrics(
            gt,
            split_path,
            crop_out / "metrics",
            f"{scan}_{crop_id}_BA_TOS_local_split",
        )

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(OUT / "local_split_manifest.csv", index=False)

    metric_rows = []
    for r in manifest_rows:
        scan = r["scan"]
        crop_id = r["crop_id"]
        mfile = OUT / scan / crop_id / "metrics" / f"{scan}_{crop_id}_BA_TOS_local_split_instance_metrics.csv"
        if not mfile.exists():
            print("[missing metric]", mfile)
            continue
        df = pd.read_csv(mfile)
        row = df.iloc[0].to_dict()
        row["scan"] = scan
        row["crop_id"] = crop_id
        metric_rows.append(row)

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "local_split_metrics_per_crop.csv", index=False)

    summary = metrics.agg({
        "precision": "mean",
        "recall": "mean",
        "f1": "mean",
        "overseg_ratio_pred_eval_gt": "mean",
        "fp": "sum",
        "fn": "sum",
        "tp": "sum",
    }).to_frame().T

    summary.insert(0, "method", "BA-TOS local-split-large-labels")
    summary.insert(1, "n", len(metrics))
    summary.to_csv(OUT / "local_split_summary.csv", index=False)

    print()
    print("==== local split manifest ====")
    print(manifest_df.to_string(index=False))

    print()
    print("==== local split summary ====")
    print(summary.to_string(index=False))

    print()
    print("wrote:", OUT / "local_split_manifest.csv")
    print("wrote:", OUT / "local_split_metrics_per_crop.csv")
    print("wrote:", OUT / "local_split_summary.csv")


if __name__ == "__main__":
    main()
