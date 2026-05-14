#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import tifffile as tiff

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(".")
DETAIL = ROOT / "results_batos/paper_tables/grain_centered_analysis/grain_centered_detail_per_gt.csv"
MANIFEST = ROOT / "results_batos/paper_tables/final_revision_clean/manual_clean_manifest_multi_iou.csv"
OUT = ROOT / "results_batos/paper_figures/grain_outlier_qc"
OUT.mkdir(parents=True, exist_ok=True)

METHOD = "BA-TOS"
N_WORST_IOU = 25
N_WORST_CENTROID = 15
HALF_SIZE = 70  # crop visual 2D around GT centroid


def read_vol(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise RuntimeError(f"expected 3D volume, got {arr.shape}: {path}")
    return arr


def find_gray_from_gt(gt_path):
    gt_path = Path(gt_path)
    parent = gt_path.parent
    candidates = list(parent.glob("*gray_crop_8u.tif")) + list(parent.glob("*gray_crop.tif"))
    if candidates:
        return candidates[0]
    return None


def centroid_of_label(labels, lab):
    coords = np.argwhere(labels == lab)
    if coords.size == 0:
        return None
    return coords.mean(axis=0)  # z,y,x


def crop2d(img2d, cy, cx, half):
    h, w = img2d.shape[:2]
    y0 = max(0, int(round(cy)) - half)
    y1 = min(h, int(round(cy)) + half)
    x0 = max(0, int(round(cx)) - half)
    x1 = min(w, int(round(cx)) + half)
    return img2d[y0:y1, x0:x1], (y0, y1, x0, x1)


def normalize_gray(slc):
    slc = slc.astype(float)
    lo, hi = np.percentile(slc, [1, 99])
    if hi <= lo:
        hi = slc.max() if slc.max() > lo else lo + 1
    out = np.clip((slc - lo) / (hi - lo), 0, 1)
    return out


def make_overlay(gray, gt_mask, pred_mask):
    g = normalize_gray(gray)
    rgb = np.stack([g, g, g], axis=-1)

    # GT in red, prediction in green, overlap in yellow
    gt_only = gt_mask & (~pred_mask)
    pred_only = pred_mask & (~gt_mask)
    overlap = gt_mask & pred_mask

    rgb[gt_only] = [1.0, 0.0, 0.0]
    rgb[pred_only] = [0.0, 0.9, 0.0]
    rgb[overlap] = [1.0, 1.0, 0.0]
    return rgb


def load_manifest():
    mf = pd.read_csv(MANIFEST)
    rows = {}
    for _, r in mf.iterrows():
        method = r.get("method", r.get("variant", None))
        key = (r["scan"], r["crop_id"], method)
        rows[key] = {
            "gt": r["gt"],
            "pred": r["pred"],
        }
    return rows


def main():
    detail = pd.read_csv(DETAIL)
    detail = detail[detail["method"] == METHOD].copy()

    worst_iou = detail.sort_values("best_iou", ascending=True).head(N_WORST_IOU)
    worst_centroid = detail.sort_values("centroid_error", ascending=False).head(N_WORST_CENTROID)

    selected = pd.concat([worst_iou, worst_centroid], ignore_index=True)
    selected = selected.drop_duplicates(subset=["scan", "crop_id", "gt_label", "best_pred_label"])
    selected = selected.sort_values(["best_iou", "centroid_error"], ascending=[True, False])

    selected_csv = OUT / "selected_outliers_for_visual_qc.csv"
    selected.to_csv(selected_csv, index=False)

    manifest = load_manifest()

    print("selected outliers:", len(selected))
    print("wrote:", selected_csv)

    # cache volumes by scan/crop
    cache = {}

    for rank, (_, r) in enumerate(selected.iterrows(), start=1):
        scan = r["scan"]
        crop_id = r["crop_id"]
        gt_label = int(r["gt_label"])
        pred_label = int(r["best_pred_label"])

        key = (scan, crop_id, METHOD)
        if key not in manifest:
            print("[skip] missing manifest", key)
            continue

        gt_path = Path(manifest[key]["gt"])
        pred_path = Path(manifest[key]["pred"])
        gray_path = find_gray_from_gt(gt_path)

        cache_key = (str(gt_path), str(pred_path), str(gray_path))
        if cache_key not in cache:
            gt = read_vol(gt_path)
            pred = read_vol(pred_path)
            gray = read_vol(gray_path) if gray_path is not None and Path(gray_path).exists() else np.zeros_like(gt, dtype=np.uint8)
            cache[cache_key] = (gray, gt, pred)
        else:
            gray, gt, pred = cache[cache_key]

        c = centroid_of_label(gt, gt_label)
        if c is None:
            print("[skip] GT label not found", gt_label)
            continue

        z, y, x = c
        z = int(round(z))
        z = max(0, min(z, gt.shape[0] - 1))

        gray_slc = gray[z]
        gt_mask = gt[z] == gt_label
        pred_mask = pred[z] == pred_label

        gray_crop, box = crop2d(gray_slc, y, x, HALF_SIZE)
        y0, y1, x0, x1 = box
        gt_crop = gt_mask[y0:y1, x0:x1]
        pred_crop = pred_mask[y0:y1, x0:x1]
        overlay = make_overlay(gray_crop, gt_crop, pred_crop)

        fig, axes = plt.subplots(1, 4, figsize=(14, 4))

        axes[0].imshow(normalize_gray(gray_crop), cmap="gray")
        axes[0].set_title("gray")

        axes[1].imshow(gt_crop, cmap="gray")
        axes[1].set_title(f"GT {gt_label}")

        axes[2].imshow(pred_crop, cmap="gray")
        axes[2].set_title(f"Pred {pred_label}")

        axes[3].imshow(overlay)
        axes[3].set_title("red=GT, green=pred, yellow=overlap")

        for ax in axes:
            ax.axis("off")

        title = (
            f"{scan} {crop_id} | GT={gt_label} Pred={pred_label} | "
            f"IoU={r['best_iou']:.3f} | "
            f"vol ratio={r['volume_ratio']:.2f} | "
            f"cent err={r['centroid_error']:.2f}"
        )
        fig.suptitle(title, fontsize=10)
        fig.tight_layout()

        out_png = OUT / f"{rank:03d}_{scan}_{crop_id}_gt{gt_label}_pred{pred_label}_iou{r['best_iou']:.3f}.png"
        fig.savefig(out_png, dpi=180)
        plt.close(fig)

        print("wrote:", out_png)

    print()
    print("done. output dir:", OUT)


if __name__ == "__main__":
    main()
