#!/usr/bin/env python3
from pathlib import Path
import csv
import numpy as np
import pandas as pd
import tifffile as tiff
from scipy import ndimage as ndi

ROOT = Path(".")
MANIFEST = ROOT / "results_batos/paper_tables/grain_centered_analysis_batos_ls/manual_manifest_batos_ls_for_grain_centered.csv"
OUT = ROOT / "results_batos/paper_tables/grain_centered_analysis_batos_ls"
OUT.mkdir(parents=True, exist_ok=True)

MIN_GT_AREA = 1000
IGNORE_GT_BORDER = True

METHOD_ORDER = ["BA-TOS-LS"]


def read_volume(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D volume, got {arr.shape}: {path}")
    return arr.astype(np.int64, copy=False)


def border_ids(labels):
    ids = set()
    faces = [
        labels[0, :, :], labels[-1, :, :],
        labels[:, 0, :], labels[:, -1, :],
        labels[:, :, 0], labels[:, :, -1],
    ]
    for face in faces:
        vals = np.unique(face)
        ids.update(int(v) for v in vals if v > 0)
    return ids


def label_stats(labels, min_area=0, ignore_border=False):
    ids, counts = np.unique(labels, return_counts=True)
    b_ids = border_ids(labels) if ignore_border else set()

    valid = []
    volume = {}

    for lab, count in zip(ids, counts):
        lab = int(lab)
        count = int(count)
        if lab == 0:
            continue
        if count < min_area:
            continue
        if lab in b_ids:
            continue
        valid.append(lab)
        volume[lab] = count

    if not valid:
        return [], volume, {}

    coms = ndi.center_of_mass(np.ones(labels.shape, dtype=np.uint8), labels, valid)
    centroids = {
        int(lab): np.array(com, dtype=float)
        for lab, com in zip(valid, coms)
    }

    return valid, volume, centroids


def all_label_volumes(labels):
    ids, counts = np.unique(labels, return_counts=True)
    return {int(i): int(c) for i, c in zip(ids, counts) if int(i) > 0}


def all_centroids(labels, ids):
    if not ids:
        return {}
    coms = ndi.center_of_mass(np.ones(labels.shape, dtype=np.uint8), labels, ids)
    return {
        int(lab): np.array(com, dtype=float)
        for lab, com in zip(ids, coms)
    }


def compute_best_matches(gt, pred, valid_gt_ids, gt_vol, gt_centroids):
    pred_vol = all_label_volumes(pred)
    pred_ids = sorted(pred_vol.keys())
    pred_centroids = all_centroids(pred, pred_ids)

    # Overlap counts between GT and predicted labels
    mask = (gt > 0) & (pred > 0)
    gt_flat = gt[mask].ravel()
    pred_flat = pred[mask].ravel()

    overlap = {}
    if len(gt_flat) > 0:
        pairs = np.stack([gt_flat, pred_flat], axis=1)
        unique_pairs, counts = np.unique(pairs, axis=0, return_counts=True)
        for (g, p), c in zip(unique_pairs, counts):
            g = int(g)
            p = int(p)
            if g in gt_vol and p > 0:
                overlap[(g, p)] = int(c)

    rows = []

    for g in valid_gt_ids:
        best_p = 0
        best_iou = 0.0
        best_inter = 0

        for (gg, pp), inter in overlap.items():
            if gg != g:
                continue
            union = gt_vol[g] + pred_vol.get(pp, 0) - inter
            iou = inter / union if union > 0 else 0.0
            if iou > best_iou:
                best_iou = iou
                best_p = pp
                best_inter = inter

        recovered_050 = int(best_iou >= 0.50)
        recovered_075 = int(best_iou >= 0.75)

        if best_p > 0 and best_p in pred_vol:
            volume_ratio = pred_vol[best_p] / gt_vol[g] if gt_vol[g] > 0 else np.nan
            relative_volume_error = abs(pred_vol[best_p] - gt_vol[g]) / gt_vol[g] if gt_vol[g] > 0 else np.nan

            if g in gt_centroids and best_p in pred_centroids:
                centroid_error = float(np.linalg.norm(gt_centroids[g] - pred_centroids[best_p]))
            else:
                centroid_error = np.nan
        else:
            volume_ratio = np.nan
            relative_volume_error = np.nan
            centroid_error = np.nan

        rows.append({
            "gt_label": int(g),
            "best_pred_label": int(best_p),
            "best_iou": float(best_iou),
            "intersection": int(best_inter),
            "gt_volume": int(gt_vol[g]),
            "pred_volume": int(pred_vol[best_p]) if best_p in pred_vol else 0,
            "volume_ratio": volume_ratio,
            "relative_volume_error": relative_volume_error,
            "centroid_error": centroid_error,
            "recovered_iou050": recovered_050,
            "recovered_iou075": recovered_075,
        })

    return pd.DataFrame(rows)


def load_manifest():
    if not MANIFEST.exists():
        raise SystemExit(f"Missing manifest: {MANIFEST}")

    rows = list(csv.DictReader(open(MANIFEST)))

    clean = []
    for r in rows:
        # expected columns: scan,crop_id,method,gt,pred,gt_found,pred_found
        if r.get("gt_found", "1") == "0":
            continue
        if r.get("pred_found", "1") == "0":
            continue

        method = r.get("method") or r.get("variant")
        if method not in METHOD_ORDER:
            continue

        clean.append({
            "scan": r["scan"],
            "crop_id": r["crop_id"],
            "method": method,
            "gt": r["gt"],
            "pred": r["pred"],
        })

    return clean


def fmt(x):
    if pd.isna(x):
        return "--"
    return f"{float(x):.3f}"


def main():
    manifest = load_manifest()

    print("manifest rows:", len(manifest))

    all_rows = []

    for item in manifest:
        scan = item["scan"]
        crop_id = item["crop_id"]
        method = item["method"]
        gt_path = Path(item["gt"])
        pred_path = Path(item["pred"])

        if not gt_path.exists():
            print("[skip missing GT]", gt_path)
            continue
        if not pred_path.exists():
            print("[skip missing PRED]", pred_path)
            continue

        print(f"[run] {scan} {crop_id} {method}")

        gt = read_volume(gt_path)
        pred = read_volume(pred_path)

        valid_gt_ids, gt_vol, gt_centroids = label_stats(
            gt,
            min_area=MIN_GT_AREA,
            ignore_border=IGNORE_GT_BORDER,
        )

        df = compute_best_matches(gt, pred, valid_gt_ids, gt_vol, gt_centroids)

        df.insert(0, "scan", scan)
        df.insert(1, "crop_id", crop_id)
        df.insert(2, "method", method)
        df.insert(3, "n_valid_gt", len(valid_gt_ids))

        all_rows.append(df)

    if not all_rows:
        raise SystemExit("No rows generated.")

    detail = pd.concat(all_rows, ignore_index=True)

    detail_csv = OUT / "grain_centered_detail_per_gt.csv"
    detail.to_csv(detail_csv, index=False)

    # Per-crop summary
    crop_summary = (
        detail.groupby(["scan", "crop_id", "method"], as_index=False)
        .agg(
            n_valid_gt=("gt_label", "count"),
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

    crop_csv = OUT / "grain_centered_summary_per_crop.csv"
    crop_summary.to_csv(crop_csv, index=False)

    # Method summary
    method_summary = (
        detail.groupby("method", as_index=False)
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

    order = {m: i for i, m in enumerate(METHOD_ORDER)}
    method_summary["order"] = method_summary["method"].map(order)
    method_summary = method_summary.sort_values("order").drop(columns=["order"])

    method_csv = OUT / "grain_centered_summary_by_method.csv"
    method_summary.to_csv(method_csv, index=False)

    # Latex table
    latex = []
    latex.append(r"\begin{table}[t]")
    latex.append(r"\centering")
    latex.append(r"\caption{Grain-centered analysis over valid reference grains.}")
    latex.append(r"\label{tab:grain_centered}")
    latex.append(r"\resizebox{\columnwidth}{!}{%")
    latex.append(r"\begin{tabular}{lrrrrr}")
    latex.append(r"\hline")
    latex.append(r"Method & Grains & Rec@.50 & Rec@.75 & Best IoU & Centroid err. \\")
    latex.append(r"\hline")

    for _, r in method_summary.iterrows():
        latex.append(
            f"{r['method']} & "
            f"{int(r['n_grains'])} & "
            f"{fmt(r['recall_iou050'])} & "
            f"{fmt(r['recall_iou075'])} & "
            f"{fmt(r['median_best_iou'])} & "
            f"{fmt(r['median_centroid_error'])} \\\\"
        )

    latex.append(r"\hline")
    latex.append(r"\end{tabular}%")
    latex.append(r"}")
    latex.append(r"\end{table}")

    latex_path = OUT / "table_grain_centered_latex.tex"
    latex_path.write_text("\n".join(latex), encoding="utf-8")

    print()
    print("wrote:", detail_csv)
    print("wrote:", crop_csv)
    print("wrote:", method_csv)
    print("wrote:", latex_path)

    print()
    print("=== GRAIN-CENTERED SUMMARY BY METHOD ===")
    print(method_summary.to_string(index=False))


if __name__ == "__main__":
    main()
