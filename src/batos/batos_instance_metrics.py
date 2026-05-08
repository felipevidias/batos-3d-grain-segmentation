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
    ap.add_argument("--min-gt-area", type=int, default=1000)
    ap.add_argument("--ignore-gt-border", action="store_true")
    ap.add_argument("--thresholds", default="0.25,0.50")
    ap.add_argument("--name", default="method")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gt = read3d(args.gt)
    pred = read3d(args.pred)

    if gt.shape != pred.shape:
        raise RuntimeError(f"Shape mismatch: gt={gt.shape}, pred={pred.shape}")

    gt_sizes_all = label_sizes(gt)
    pred_sizes_all = label_sizes(pred)

    gt_border = border_ids(gt) if args.ignore_gt_border else set()

    valid_gt = {
        gid for gid, sz in gt_sizes_all.items()
        if sz >= args.min_gt_area and gid not in gt_border
    }

    ignored_gt = set(gt_sizes_all.keys()) - valid_gt

    # Intersections between GT and prediction.
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

    # Decide which predicted labels are evaluable.
    # Predictions mostly overlapping ignored border GT are ignored, not counted as FP.
    pred_eval = set()
    pred_ignored = set()

    for pid, psz in pred_sizes_all.items():
        overlaps = [(gid, inter) for (gid, p), inter in intersections.items() if p == pid]

        if not overlaps:
            # pure background prediction counts as evaluable false positive
            pred_eval.add(pid)
            continue

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
        gsz = gt_sizes_all[gid]
        psz = pred_sizes_all[pid]
        union = gsz + psz - inter
        if union > 0:
            iou[gt_index[gid], pred_index[pid]] = inter / union

    thresholds = [float(x) for x in args.thresholds.split(",")]

    rows = []
    match_rows = []

    for thr in thresholds:
        if iou.size == 0:
            matched = []
        else:
            # Maximize IoU via Hungarian assignment.
            cost = 1.0 - iou
            rr, cc = linear_sum_assignment(cost)
            matched = [
                (gt_list[r], pred_list[c], float(iou[r, c]))
                for r, c in zip(rr, cc)
                if iou[r, c] >= thr
            ]

        matched_gt = {m[0] for m in matched}
        matched_pred = {m[1] for m in matched}

        tp = len(matched)
        fp = len(pred_eval) - len(matched_pred)
        fn = len(valid_gt) - len(matched_gt)

        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

        rows.append({
            "method": args.name,
            "threshold": thr,
            "gt_valid": len(valid_gt),
            "pred_total": len(pred_sizes_all),
            "pred_eval": len(pred_eval),
            "pred_ignored_border": len(pred_ignored),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "overseg_ratio_pred_eval_gt": len(pred_eval) / len(valid_gt) if valid_gt else 0.0,
        })

        for gid, pid, val in matched:
            match_rows.append({
                "threshold": thr,
                "gt_id": gid,
                "pred_id": pid,
                "iou": val,
            })

    out_csv = out_dir / f"{args.name}_instance_metrics.csv"
    with open(out_csv, "w", newline="") as f:
        fieldnames = [
            "method",
            "threshold",
            "gt_valid",
            "pred_total",
            "pred_eval",
            "pred_ignored_border",
            "tp",
            "fp",
            "fn",
            "precision",
            "recall",
            "f1",
            "overseg_ratio_pred_eval_gt",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    out_matches = out_dir / f"{args.name}_matches.csv"
    with open(out_matches, "w", newline="") as f:
        fieldnames = ["threshold", "gt_id", "pred_id", "iou"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(match_rows)

    print("wrote:", out_csv)
    print("wrote:", out_matches)
    print("")
    for r in rows:
        print(
            f"{args.name} @ IoU {r['threshold']:.2f}: "
            f"GT={r['gt_valid']}, pred_total={r['pred_total']}, pred_eval={r['pred_eval']}, "
            f"TP={r['tp']}, FP={r['fp']}, FN={r['fn']}, "
            f"P={r['precision']:.4f}, R={r['recall']:.4f}, F1={r['f1']:.4f}, "
            f"overseg={r['overseg_ratio_pred_eval_gt']:.3f}"
        )


if __name__ == "__main__":
    main()
