"""Stage 5: object-level evaluation (seed-ignore, paper §III-C) + pooled / crop-macro summaries.

Valid reference = GT label with >= 1000 voxels that touches none of the 6 crop faces; every other
GT label is "excluded".

seed_ignore (primary, reconstructed e09b, BG11):
  a prediction is IGNORED iff all of its seed (marker) voxels lie inside excluded reference
  components; a prediction seeded in background or in a valid reference is evaluated (FP if unmatched).
  IoU is computed on the domain without excluded-reference voxels (eval.domain=exclude_ignored_gt,
  as e12) or on the full crop (eval.domain=full).
strict (secondary, exact port of e12_boundary_margin_audit at margin 0):
  every prediction with >= 1 voxel outside excluded references is evaluated.
Matching: maximum-cardinality bipartite matching over pairs with IoU >= threshold (as e12).

Pooled-object and crop-macro summaries are written to separate, labelled files (CLAUDE.md rule 5).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

from .common import load_config, read3d


def reference_sets(gt: np.ndarray, min_voxels: int, exclude_face_touching: bool = True):
    gt = gt.astype(np.int64, copy=False)
    ids, counts = np.unique(gt, return_counts=True)
    keep = ids > 0
    ids, counts = ids[keep], counts[keep]
    faces = np.concatenate([gt[0].ravel(), gt[-1].ravel(), gt[:, 0].ravel(), gt[:, -1].ravel(),
                            gt[:, :, 0].ravel(), gt[:, :, -1].ravel()])
    border = set(np.unique(faces[faces > 0]).tolist()) if exclude_face_touching else set()
    valid = {int(i) for i, c in zip(ids, counts) if c >= min_voxels and int(i) not in border}
    excluded = {int(i) for i in ids} - valid
    return valid, excluded


def _match(pred_ids, gt_ids, pairs_iou, threshold):
    """pairs_iou: dict (pred_idx, gt_idx) -> iou. Returns matched IoU list."""
    selected = [(p, g, v) for (p, g), v in pairs_iou.items() if v >= threshold]
    if not selected or not len(pred_ids) or not len(gt_ids):
        return []
    rows = np.fromiter((s[0] for s in selected), dtype=np.int32)
    cols = np.fromiter((s[1] for s in selected), dtype=np.int32)
    graph = csr_matrix((np.ones(len(selected), dtype=np.uint8), (rows, cols)), shape=(len(pred_ids), len(gt_ids)))
    row_to_col = maximum_bipartite_matching(graph, perm_type="column")
    lookup = {(p, g): v for p, g, v in selected}
    return [lookup[(r, int(c))] for r, c in enumerate(row_to_col) if c >= 0 and (r, int(c)) in lookup]


def evaluate_protocol(pred: np.ndarray, gt: np.ndarray, valid: set[int], excluded: set[int],
                      eval_pred_ids: list[int], domain: np.ndarray, thresholds) -> dict[str, Any]:
    gt_ids = sorted(valid)
    pred_ids = sorted(eval_pred_ids)
    gt_index = {g: i for i, g in enumerate(gt_ids)}
    pred_index = {p: i for i, p in enumerate(pred_ids)}
    pred_d = pred[domain].astype(np.int64)
    gt_d = gt[domain].astype(np.int64)
    pred_area = pd.Series(pred_d[pred_d > 0]).value_counts()
    gt_area = pd.Series(gt_d[gt_d > 0]).value_counts()
    both = (pred_d > 0) & (gt_d > 0)
    pairs = pd.DataFrame({"p": pred_d[both], "g": gt_d[both]}).value_counts()
    pairs_iou = {}
    for (p, g), inter in pairs.items():
        if p in pred_index and g in gt_index:
            union = int(pred_area[p]) + int(gt_area[g]) - int(inter)
            pairs_iou[(pred_index[p], gt_index[g])] = int(inter) / max(union, 1)
    out: dict[str, Any] = {"n_gt_valid": len(gt_ids), "n_pred_eval": len(pred_ids)}
    for thr in thresholds:
        suffix = f"{int(round(thr * 100)):03d}"
        matched = _match(pred_ids, gt_ids, pairs_iou, thr)
        tp = len(matched)
        fp, fn = len(pred_ids) - tp, len(gt_ids) - tp
        precision = tp / len(pred_ids) if pred_ids else 0.0
        recall = tp / len(gt_ids) if gt_ids else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        out.update({f"tp_iou{suffix}": tp, f"fp_iou{suffix}": fp, f"fn_iou{suffix}": fn,
                    f"precision_iou{suffix}": precision, f"recall_iou{suffix}": recall, f"f1_iou{suffix}": f1,
                    f"mean_matched_iou{suffix}": float(np.mean(matched)) if matched else 0.0})
    out["overseg_ratio"] = len(pred_ids) / len(gt_ids) if gt_ids else 0.0
    return out


def evaluate_crop(pred: np.ndarray, seeds: np.ndarray, gt: np.ndarray, cfg: dict[str, Any]) -> list[dict[str, Any]]:
    if not (pred.shape == seeds.shape == gt.shape):
        raise ValueError(f"shape mismatch pred={pred.shape} seeds={seeds.shape} gt={gt.shape}")
    pred = pred.astype(np.int64, copy=False)
    gt = gt.astype(np.int64, copy=False)
    valid, excluded = reference_sets(gt, cfg["min_gt_voxels"], cfg["exclude_face_touching"])
    excluded_mask = np.isin(gt, list(excluded)) if excluded else np.zeros(gt.shape, bool)
    pred_ids = [int(p) for p in np.unique(pred) if p > 0]
    thresholds = cfg["iou_thresholds"]

    # seed_ignore: seeds of prediction p = seed voxels whose basin label is p
    seed_vox = (seeds > 0) & (pred > 0)
    seed_pred = pred[seed_vox]
    seed_in_excluded = excluded_mask[seed_vox]
    frame = pd.DataFrame({"p": seed_pred, "ex": seed_in_excluded})
    all_excluded = frame.groupby("p")["ex"].all() if len(frame) else pd.Series(dtype=bool)
    ignored = {int(p) for p, flag in all_excluded.items() if flag}
    unseeded = sorted(set(pred_ids) - set(all_excluded.index.astype(int)))
    eval_ids = [p for p in pred_ids if p not in ignored]
    domain = ~excluded_mask if cfg["domain"] == "exclude_ignored_gt" else np.ones(gt.shape, bool)
    seed_row = {"evaluation_protocol": "seed_ignore", "n_pred_total": len(pred_ids),
                "n_pred_ignored": len(ignored), "n_pred_unseeded": len(unseeded),
                **evaluate_protocol(pred, gt, valid, excluded, eval_ids, domain, thresholds)}

    strict_ids = [int(p) for p in np.unique(pred[~excluded_mask]) if p > 0]
    strict_row = {"evaluation_protocol": "strict", "n_pred_total": len(pred_ids),
                  "n_pred_ignored": len(pred_ids) - len(strict_ids), "n_pred_unseeded": len(unseeded),
                  **evaluate_protocol(pred, gt, valid, excluded, strict_ids, ~excluded_mask, thresholds)}
    return [seed_row, strict_row]


def summarize(per_crop: pd.DataFrame, thresholds) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (pooled_object, crop_macro) summaries grouped by method / size_filter / protocol."""
    keys = ["method", "size_filter", "evaluation_protocol"]
    pooled_rows, macro_rows = [], []
    for key, group in per_crop.groupby(keys, sort=False):
        pooled = dict(zip(keys, key), summary_type="pooled_object", n_crops=len(group),
                      n_gt_valid=int(group["n_gt_valid"].sum()), n_pred_eval=int(group["n_pred_eval"].sum()))
        macro = dict(zip(keys, key), summary_type="crop_macro", n_crops=len(group))
        for thr in thresholds:
            s = f"{int(round(thr * 100)):03d}"
            tp, fp, fn = (int(group[f"{m}_iou{s}"].sum()) for m in ("tp", "fp", "fn"))
            p = tp / (tp + fp) if tp + fp else 0.0
            r = tp / (tp + fn) if tp + fn else 0.0
            pooled.update({f"tp_iou{s}": tp, f"fp_iou{s}": fp, f"fn_iou{s}": fn, f"precision_iou{s}": p,
                           f"recall_iou{s}": r, f"f1_iou{s}": 2 * p * r / (p + r) if p + r else 0.0})
            f1 = group[f"f1_iou{s}"].astype(float)
            macro.update({f"f1_iou{s}_mean": f1.mean(), f"f1_iou{s}_sd_ddof1": f1.std(ddof=1) if len(f1) > 1 else 0.0,
                          f"f1_iou{s}_sd_ddof0": f1.std(ddof=0), f"f1_iou{s}_median": f1.median(),
                          f"f1_iou{s}_min": f1.min()})
        pooled_rows.append(pooled)
        macro_rows.append(macro)
    return pd.DataFrame(pooled_rows), pd.DataFrame(macro_rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="batos-paper-v2 stage 5: evaluation")
    sub = parser.add_subparsers(dest="command", required=True)
    one = sub.add_parser("crop", help="evaluate one segmentation")
    one.add_argument("--pred", type=Path, required=True)
    one.add_argument("--seeds", type=Path, required=True, help="normalized markers used by the watershed")
    one.add_argument("--gt", type=Path, required=True)
    one.add_argument("--crop-id", required=True)
    one.add_argument("--method", required=True)
    one.add_argument("--size-filter", default="off")
    one.add_argument("--out", type=Path, required=True)
    agg = sub.add_parser("aggregate", help="pooled-object and crop-macro summaries")
    agg.add_argument("--inputs", type=Path, nargs="+", required=True, help="per-crop CSVs")
    agg.add_argument("--out-dir", type=Path, required=True)
    for p in (one, agg):
        p.add_argument("--config", type=Path, default=Path("configs/batos_paper_v2.yaml"))
    args = parser.parse_args(argv)
    cfg = load_config(args.config)["eval"]
    if args.command == "crop":
        rows = evaluate_crop(read3d(args.pred), read3d(args.seeds), read3d(args.gt), cfg)
        frame = pd.DataFrame([{"crop_id": args.crop_id, "method": args.method, "size_filter": args.size_filter, **r}
                              for r in rows])
        args.out.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.out, index=False)
        print(frame[["evaluation_protocol", "f1_iou025", "f1_iou050", "f1_iou075"]].to_string(index=False))
    else:
        per_crop = pd.concat([pd.read_csv(p) for p in args.inputs], ignore_index=True)
        pooled, macro = summarize(per_crop, cfg["iou_thresholds"])
        args.out_dir.mkdir(parents=True, exist_ok=True)
        per_crop.to_csv(args.out_dir / "per_crop_metrics.csv", index=False)
        pooled.to_csv(args.out_dir / "summary_pooled_object.csv", index=False)
        macro.to_csv(args.out_dir / "summary_crop_macro.csv", index=False)
        print(f"wrote {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
