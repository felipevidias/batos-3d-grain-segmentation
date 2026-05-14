#!/usr/bin/env python3
from pathlib import Path
import csv
import subprocess
import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi
import pandas as pd

ROOT = Path(".")
MANIFEST = ROOT / "results_batos/paper_tables/final_revision_clean/manual_clean_manifest_multi_iou.csv"
OUT = ROOT / "results_batos/paper_tables/batos_split_disconnected"
OUT.mkdir(parents=True, exist_ok=True)

PY = "/home/felipe/Pesquisa-Grain_Seg/Filtro_Conexo/venv/bin/python3"
METRICS = ROOT / "code/segmentation/batos/batos_instance_metrics.py"

METHOD = "BA-TOS"
MIN_COMPONENT_SIZE = 500
CONNECTIVITY = 1


def read_labels(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D label volume, got {arr.shape}: {path}")
    return arr.astype(np.int32, copy=False)


def split_disconnected(labels, min_size=500):
    structure = ndi.generate_binary_structure(3, CONNECTIVITY)
    out = np.zeros_like(labels, dtype=np.int32)

    next_id = 1
    split_count = 0
    removed_count = 0
    original_count = 0
    new_count = 0

    ids = np.unique(labels)
    ids = ids[ids > 0]

    for lab in ids:
        original_count += 1
        mask = labels == lab
        cc, ncc = ndi.label(mask, structure=structure)

        kept_for_this = 0
        for cid in range(1, ncc + 1):
            comp = cc == cid
            size = int(comp.sum())

            if size < min_size:
                removed_count += 1
                continue

            out[comp] = next_id
            next_id += 1
            kept_for_this += 1
            new_count += 1

        if kept_for_this > 1:
            split_count += 1

    return out, {
        "original_labels": int(original_count),
        "new_labels": int(new_count),
        "labels_split": int(split_count),
        "small_components_removed": int(removed_count),
    }


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
    rows = load_manifest()
    print("BA-TOS rows:", len(rows))

    split_rows = []

    for r in rows:
        scan = r["scan"]
        crop_id = r["crop_id"]
        gt = Path(r["gt"])
        pred = Path(r["pred"])

        labels = read_labels(pred)
        split_labels, stats = split_disconnected(labels, MIN_COMPONENT_SIZE)

        crop_out = OUT / scan / crop_id
        crop_out.mkdir(parents=True, exist_ok=True)

        split_path = crop_out / f"{scan}_{crop_id}_batos_split_disconnected_labels.tif"
        tiff.imwrite(str(split_path), split_labels.astype(np.int32))

        stats_row = {
            "scan": scan,
            "crop_id": crop_id,
            "gt": str(gt),
            "original_pred": str(pred),
            "split_pred": str(split_path),
            **stats,
        }
        split_rows.append(stats_row)

        run_metrics(
            gt,
            split_path,
            crop_out / "metrics",
            f"{scan}_{crop_id}_BA_TOS_split_disconnected",
        )

    stats_df = pd.DataFrame(split_rows)
    stats_df.to_csv(OUT / "split_disconnected_manifest.csv", index=False)

    # collect metrics
    metric_rows = []
    for r in split_rows:
        scan = r["scan"]
        crop_id = r["crop_id"]
        mfile = OUT / scan / crop_id / "metrics" / f"{scan}_{crop_id}_BA_TOS_split_disconnected_instance_metrics.csv"
        if not mfile.exists():
            print("[missing metric]", mfile)
            continue
        df = pd.read_csv(mfile)
        row = df.iloc[0].to_dict()
        row["scan"] = scan
        row["crop_id"] = crop_id
        metric_rows.append(row)

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "split_disconnected_metrics_per_crop.csv", index=False)

    summary = metrics.agg({
        "precision": "mean",
        "recall": "mean",
        "f1": "mean",
        "overseg_ratio_pred_eval_gt": "mean",
        "fp": "sum",
        "fn": "sum",
        "tp": "sum",
    }).to_frame().T

    summary.insert(0, "method", "BA-TOS split-disconnected")
    summary.insert(1, "n", len(metrics))
    summary.to_csv(OUT / "split_disconnected_summary.csv", index=False)

    print()
    print("==== split stats ====")
    print(stats_df.to_string(index=False))

    print()
    print("==== summary ====")
    print(summary.to_string(index=False))

    print()
    print("wrote:", OUT / "split_disconnected_manifest.csv")
    print("wrote:", OUT / "split_disconnected_metrics_per_crop.csv")
    print("wrote:", OUT / "split_disconnected_summary.csv")


if __name__ == "__main__":
    main()
