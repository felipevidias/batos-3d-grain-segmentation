#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
import numpy as np
import tifffile as tiff
from scipy.optimize import linear_sum_assignment


def read3d(path):
    arr = np.asarray(tiff.imread(str(path)))
    if arr.ndim == 4 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D volume, got {arr.shape} for {path}")
    return arr.astype(np.int64)


def border_ids(vol):
    vals = np.concatenate([
        vol[0, :, :].ravel(),
        vol[-1, :, :].ravel(),
        vol[:, 0, :].ravel(),
        vol[:, -1, :].ravel(),
        vol[:, :, 0].ravel(),
        vol[:, :, -1].ravel(),
    ])
    vals = vals[vals > 0]
    return set(vals.astype(np.int64).tolist())


def label_sizes(vol):
    ids, counts = np.unique(vol, return_counts=True)
    keep = ids > 0
    return {int(i): int(c) for i, c in zip(ids[keep], counts[keep])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gt", required=True)
    ap.add_argument("--pred", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--name", default="diag")
    ap.add_argument("--min-gt-area", type=int, default=1000)
    ap.add_argument("--threshold", type=float, default=0.50)
    ap.add_argument("--ignore-gt-border", action="store_true")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gt = read3d(args.gt)
    pred = read3d(args.pred)

    if gt.shape != pred.shape:
        raise RuntimeError(f"Shape mismatch: gt={gt.shape}, pred={pred.shape}")

    gt_sizes = label_sizes(gt)
    pred_sizes = label_sizes(pred)

    gt_border = border_ids(gt) if args.ignore_gt_border else set()

    valid_gt = {
        gid for gid, sz in gt_sizes.items()
        if sz >= args.min_gt_area and gid not in gt_border
    }

    ignored_gt = set(gt_sizes.keys()) - valid_gt

    pair_mask = (gt > 0) & (pred > 0)
    gt_pair = gt[pair_mask]
    pred_pair = pred[pair_mask]

    intersections = {}
    if gt_pair.size:
        max_pred = int(pred.max()) + 1
        codes = gt_pair * max_pred + pred_pair
        codes_u, counts = np.unique(codes, return_counts=True)
        for code, cnt in zip(codes_u, counts):
            gid = int(code // max_pred)
            pid = int(code % max_pred)
            intersections[(gid, pid)] = int(cnt)

    pred_eval = set()
    pred_ignored = set()

    for pid, psz in pred_sizes.items():
        overlaps = [(gid, inter) for (gid, p), inter in intersections.items() if p == pid]
        valid_overlap = sum(inter for gid, inter in overlaps if gid in valid_gt)
        ignored_overlap = sum(inter for gid, inter in overlaps if gid in ignored_gt)

        if ignored_overlap > valid_overlap:
            pred_ignored.add(pid)
        else:
            pred_eval.add(pid)

    gt_list = sorted(valid_gt)
    pred_list = sorted(pred_eval)

    gt_index = {gid: i for i, gid in enumerate(gt_list)}
    pred_index = {pid: j for j, pid in enumerate(pred_list)}

    iou = np.zeros((len(gt_list), len(pred_list)), dtype=np.float64)

    for (gid, pid), inter in intersections.items():
        if gid not in valid_gt or pid not in pred_eval:
            continue
        union = gt_sizes[gid] + pred_sizes[pid] - inter
        if union > 0:
            iou[gt_index[gid], pred_index[pid]] = inter / union

    if iou.size:
        rr, cc = linear_sum_assignment(1.0 - iou)
        matches = [
            (gt_list[r], pred_list[c], float(iou[r, c]))
            for r, c in zip(rr, cc)
            if iou[r, c] >= args.threshold
        ]
    else:
        matches = []

    matched_gt = {m[0] for m in matches}
    matched_pred = {m[1] for m in matches}

    fp_pred = sorted(pred_eval - matched_pred)
    fn_gt = sorted(valid_gt - matched_gt)

    fp_rows = []

    for pid in fp_pred:
        overlaps = [(gid, inter) for (gid, p), inter in intersections.items() if p == pid]
        pred_size = pred_sizes[pid]

        valid_overlap = sum(inter for gid, inter in overlaps if gid in valid_gt)
        ignored_overlap = sum(inter for gid, inter in overlaps if gid in ignored_gt)
        bg_overlap = pred_size - valid_overlap - ignored_overlap

        best_gid = 0
        best_iou = 0.0
        best_inter = 0

        for gid, inter in overlaps:
            if gid not in valid_gt:
                continue
            union = gt_sizes[gid] + pred_size - inter
            val = inter / union if union > 0 else 0.0
            if val > best_iou:
                best_iou = val
                best_gid = gid
                best_inter = inter

        if valid_overlap == 0 and ignored_overlap > 0:
            reason = "mostly_ignored_or_border_gt"
        elif valid_overlap == 0 and bg_overlap > 0:
            reason = "background_false_positive"
        elif best_iou < args.threshold:
            reason = "fragment_or_low_iou"
        else:
            reason = "unmatched_assignment"

        fp_rows.append({
            "pred_id": pid,
            "pred_voxels": pred_size,
            "best_valid_gt": best_gid,
            "best_iou": best_iou,
            "best_intersection": best_inter,
            "valid_gt_overlap": valid_overlap,
            "ignored_gt_overlap": ignored_overlap,
            "background_overlap": bg_overlap,
            "reason": reason,
        })

    # Quantos preds tocam cada GT válido? Isso mede split.
    split_rows = []
    for gid in sorted(valid_gt):
        touched = []
        for (g, pid), inter in intersections.items():
            if g != gid or pid not in pred_eval:
                continue
            union = gt_sizes[gid] + pred_sizes[pid] - inter
            val = inter / union if union > 0 else 0.0
            touched.append((pid, inter, val))
        touched_sorted = sorted(touched, key=lambda x: x[2], reverse=True)
        split_rows.append({
            "gt_id": gid,
            "gt_voxels": gt_sizes[gid],
            "num_pred_touching": len(touched_sorted),
            "best_pred": touched_sorted[0][0] if touched_sorted else 0,
            "best_iou": touched_sorted[0][2] if touched_sorted else 0.0,
            "sum_intersections": sum(x[1] for x in touched_sorted),
        })

    fp_out = out_dir / f"{args.name}_false_positives.csv"
    with open(fp_out, "w", newline="") as f:
        fieldnames = [
            "pred_id",
            "pred_voxels",
            "best_valid_gt",
            "best_iou",
            "best_intersection",
            "valid_gt_overlap",
            "ignored_gt_overlap",
            "background_overlap",
            "reason",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(fp_rows)

    split_out = out_dir / f"{args.name}_split_by_gt.csv"
    with open(split_out, "w", newline="") as f:
        fieldnames = [
            "gt_id",
            "gt_voxels",
            "num_pred_touching",
            "best_pred",
            "best_iou",
            "sum_intersections",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(split_rows)

    summary_out = out_dir / f"{args.name}_error_summary.csv"

    from collections import Counter
    reason_counts = Counter(row["reason"] for row in fp_rows)

    with open(summary_out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["gt_valid", len(valid_gt)])
        w.writerow(["pred_total", len(pred_sizes)])
        w.writerow(["pred_eval", len(pred_eval)])
        w.writerow(["pred_ignored", len(pred_ignored)])
        w.writerow(["matched_gt", len(matched_gt)])
        w.writerow(["matched_pred", len(matched_pred)])
        w.writerow(["false_positive_count", len(fp_pred)])
        w.writerow(["false_negative_count", len(fn_gt)])
        for k, v in reason_counts.items():
            w.writerow([f"fp_reason_{k}", v])

    print("wrote:", summary_out)
    print("wrote:", fp_out)
    print("wrote:", split_out)
    print("")
    print("False positive reasons:")
    for k, v in reason_counts.items():
        print(k, v)


if __name__ == "__main__":
    main()
