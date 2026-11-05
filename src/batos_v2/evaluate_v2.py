"""Phase 1d: label-agnostic object evaluator (paper Sec. III-C, implemented geometrically).

Why: the frozen e09b_ignore_aware_audit.py decides seed-ignore by assuming that a prediction's label
value equals the label of the marker component that seeded it. That holds for the E11 watershed, but
not for methods without markers (Cellpose, SAM, nnInteractive) or after any relabelling. Here the link
between a prediction and its support is purely geometric, so the result is invariant to label values.

Reference protocol (identical to e09b):
  valid GT    = GT label with >= min_gt_voxels voxels that touches none of the 6 crop faces;
  excluded GT = every other GT label; ignore_mask = voxels of excluded GT.
Support region S(p) of prediction p:
  with markers : the marker voxels inside p ("seed region"; markers > 0 & pred == p);
  without      : p's own voxels.
Ignore rule (prediction p is removed from the evaluation, never from the segmentation):
  dominant : the most frequent GT label among the GT-positive voxels of S(p) (ties -> smallest label)
             is an excluded label. This is e09b's per-marker rule, the default with markers.
  fraction : |S(p) & ignore_mask| / |S(p)| >= --ignore-overlap-threshold. Default without markers; with
             S(p) = own voxels and threshold 0.50 it equals e09b's "overlap_ignore" protocol.
  A support with no GT-positive voxel ("background-seeded") is never ignored: unmatched -> FP.
  With markers, a prediction containing no marker voxel ("unseeded") is evaluated (as in e09b).
Metrics (on the domain ~ignore_mask, ignored predictions zeroed, as e09b):
  P / R / F1 / TP / FP / FN / mean matched IoU at IoU 0.25, 0.50, 0.75. Matching = maximum-cardinality
  bipartite matching over pairs with IoU >= t (TP identical to e09b); among maximum matchings the one with
  maximum total IoU is taken, so mean matched IoU is label-invariant. e09b keeps scipy's arbitrary maximum
  matching, so its mean_matched_iou025 depends on label order; at t >= 0.50 the matching is unique and
  everything equals e09b;
  ap_iouXXX      = TP / (TP + FP + FN) at IoU t (Cellpose / StarDist "average precision"; no scores);
  coco_ap_iouXXX = score-ranked 101-point interpolated AP (greedy matching by score), only when
                   per-prediction scores are supplied, else NaN;
  morphology of objects matched at --morph-iou (default 0.50), computed on the full objects (not the
  evaluation domain): volume, equivalent diameter, centroid, surface area (marching cubes), sphericity,
  elongation (major / minor principal axis length), and their errors (pred - GT).
Pooled-object and crop-macro summaries are NOT mixed here; this module produces per-crop rows only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.sparse import csr_matrix
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import connected_components

THRESHOLDS = (0.25, 0.50, 0.75)


def _suffix(t: float) -> str:
    return f"{int(round(t * 100)):03d}"


# --------------------------------------------------------------------------- reference protocol
def reference_protocol(gt: np.ndarray, min_gt_voxels: int = 1000, exclude_face_touching: bool = True) -> dict:
    gt = np.asarray(gt)
    if gt.ndim != 3:
        raise ValueError(f"GT must be 3D, got {gt.shape}")
    if (gt < 0).any():
        raise ValueError("GT contains negative labels")
    gt = gt.astype(np.int64, copy=False)
    areas = np.bincount(gt.ravel())
    border = set()
    if exclude_face_touching:
        faces = np.concatenate([gt[0].ravel(), gt[-1].ravel(), gt[:, 0].ravel(), gt[:, -1].ravel(),
                                gt[:, :, 0].ravel(), gt[:, :, -1].ravel()])
        border = {int(v) for v in np.unique(faces) if v > 0}
    valid_ids = np.asarray([lab for lab in range(1, len(areas))
                            if areas[lab] >= min_gt_voxels and lab not in border], dtype=np.int64)
    present = np.flatnonzero(areas)
    excluded_ids = np.asarray([lab for lab in present if lab > 0 and lab not in set(valid_ids.tolist())],
                              dtype=np.int64)
    valid_lookup = np.zeros(len(areas), dtype=bool)
    valid_lookup[valid_ids] = True
    remap = np.zeros(len(areas), dtype=np.int64)
    remap[valid_ids] = np.arange(1, len(valid_ids) + 1)
    valid_mask = valid_lookup[gt]
    return {"gt": gt, "valid_ids": valid_ids, "excluded_ids": excluded_ids, "valid_lookup": valid_lookup,
            "valid_gt": np.where(valid_mask, remap[gt], 0), "ignore_mask": (gt > 0) & ~valid_mask}


# --------------------------------------------------------------------------- ignore decision
def classify_predictions(pred: np.ndarray, protocol: dict, markers: np.ndarray | None = None,
                         rule: str | None = None, threshold: float = 0.50) -> pd.DataFrame:
    """One row per prediction label: support size, GT composition, category and the ignore decision."""
    if rule is None:
        rule = "dominant" if markers is not None else "fraction"
    if rule not in ("dominant", "fraction"):
        raise ValueError(f"unknown rule {rule!r}")
    gt = protocol["gt"]
    if pred.shape != gt.shape or (markers is not None and markers.shape != gt.shape):
        raise ValueError(f"shape mismatch pred={pred.shape} gt={gt.shape} "
                         f"markers={None if markers is None else markers.shape}")
    pred = pred.astype(np.int64, copy=False)
    pred_ids, pred_vox = np.unique(pred[pred > 0], return_counts=True)
    support = (pred > 0) if markers is None else ((pred > 0) & (markers > 0))
    sp, sg = pred[support], gt[support]
    s_ex = protocol["ignore_mask"][support]

    table = pd.DataFrame({"prediction_label": pred_ids, "prediction_voxels": pred_vox})
    n_sup = pd.Series(sp).value_counts()
    n_pos = pd.Series(sp[sg > 0]).value_counts()
    n_exc = pd.Series(sp[s_ex]).value_counts()
    table["support_voxels"] = table["prediction_label"].map(n_sup).fillna(0).astype(np.int64)
    table["support_gt_positive_voxels"] = table["prediction_label"].map(n_pos).fillna(0).astype(np.int64)
    table["support_excluded_voxels"] = table["prediction_label"].map(n_exc).fillna(0).astype(np.int64)
    table["support_excluded_fraction"] = (table["support_excluded_voxels"]
                                          / table["support_voxels"].where(table["support_voxels"] > 0))
    # dominant GT label among GT-positive support voxels; ties -> smallest label (as np.unique + argmax)
    pos = sg > 0
    if pos.any():
        pairs = pd.DataFrame({"p": sp[pos], "g": sg[pos]}).value_counts().rename("n").reset_index()
        pairs = pairs.sort_values(["p", "n", "g"], ascending=[True, False, True]).drop_duplicates("p")
        dominant = pairs.set_index("p")["g"]
    else:
        dominant = pd.Series(dtype=np.int64)
    table["dominant_gt_label"] = table["prediction_label"].map(dominant).fillna(0).astype(np.int64)
    lookup = protocol["valid_lookup"]
    dom = table["dominant_gt_label"].to_numpy()
    dom_valid = np.zeros(len(dom), dtype=bool)
    inside = (dom > 0) & (dom < len(lookup))
    dom_valid[inside] = lookup[dom[inside]]
    category = np.where(table["support_voxels"] == 0, "unseeded",
                        np.where(dom == 0, "background", np.where(dom_valid, "valid_gt", "ignored_gt")))
    table["support_category"] = category
    if rule == "dominant":
        ignored = category == "ignored_gt"
    else:
        ignored = (table["support_excluded_fraction"].fillna(0.0) >= threshold).to_numpy()
    table["ignored"] = ignored
    table["support"] = "seed_region" if markers is not None else "own_voxels"
    table["rule"] = rule
    table["threshold"] = threshold if rule == "fraction" else np.nan
    return table


# --------------------------------------------------------------------------- matching / metrics
def _iou_edges(pred: np.ndarray, protocol: dict, eval_labels: np.ndarray):
    """IoU of every (evaluated prediction, valid GT) pair on the domain ~ignore_mask (as e09b)."""
    domain = ~protocol["ignore_mask"]
    pv = pred[domain].astype(np.int64, copy=False)
    gv = protocol["valid_gt"][domain]
    keep = np.isin(pv, eval_labels)
    pv = np.where(keep, pv, 0)
    index = np.zeros(int(pred.max()) + 1 if pred.size else 1, dtype=np.int64)
    index[eval_labels] = np.arange(1, len(eval_labels) + 1)
    pc = index[pv]
    n_pred, n_gt = len(eval_labels), len(protocol["valid_ids"])
    pred_areas = np.bincount(pc, minlength=n_pred + 1)
    gt_areas = np.bincount(gv, minlength=n_gt + 1)
    both = (pc > 0) & (gv > 0)
    code = pc[both] * (n_gt + 1) + gv[both]
    pairs, inter = np.unique(code, return_counts=True)
    pp, gg = pairs // (n_gt + 1), pairs % (n_gt + 1)
    iou = inter / np.maximum(pred_areas[pp] + gt_areas[gg] - inter, 1)
    return pp - 1, gg - 1, iou, n_pred, n_gt


def _max_card_matching(pp, gg, iou, n_pred, n_gt, t):
    """Maximum-cardinality matching over pairs with IoU >= t (TP identical to e09b), and among those the
    one with maximum total IoU, so the matched set does not depend on label order. e09b keeps scipy's
    arbitrary maximum matching instead; the two can differ only where a node has >= 2 candidate edges,
    i.e. only for t < 0.5 (for t >= 0.5 each prediction / GT has at most one partner with IoU > 0.5)."""
    empty = (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0))
    sel = iou >= t
    if n_pred == 0 or n_gt == 0 or not sel.any():
        return empty
    ep, eg, ev = pp[sel], gg[sel], iou[sel]
    # connected components of the candidate graph (pred nodes 0..n_pred-1, GT nodes n_pred..)
    graph = csr_matrix((np.ones(len(ep), np.uint8), (ep, eg + n_pred)), shape=(n_pred + n_gt,) * 2)
    _, comp = connected_components(graph, directed=False)
    out_p, out_g, out_v = [], [], []
    order = np.argsort(comp[ep], kind="stable")
    bounds = np.flatnonzero(np.diff(comp[ep][order])) + 1
    for idx in np.split(order, bounds):
        if len(idx) == 1:
            out_p.append(ep[idx]), out_g.append(eg[idx]), out_v.append(ev[idx])
            continue
        up, ip = np.unique(ep[idx], return_inverse=True)
        ug, ig = np.unique(eg[idx], return_inverse=True)
        big = float(min(len(up), len(ug)) + 1)  # cardinality dominates total IoU
        w = np.zeros((len(up), len(ug)))
        w[ip, ig] = big + ev[idx]
        r, c = linear_sum_assignment(w, maximize=True)
        keep = w[r, c] > 0
        out_p.append(up[r[keep]]), out_g.append(ug[c[keep]]), out_v.append(w[r[keep], c[keep]] - big)
    return np.concatenate(out_p), np.concatenate(out_g), np.concatenate(out_v)


def coco_ap(pp, gg, iou, scores: np.ndarray, n_gt: int, t: float) -> float:
    """101-point interpolated AP: predictions ranked by score (desc, stable), each greedily takes the
    unmatched valid GT of highest IoU >= t."""
    if n_gt == 0:
        return float("nan")
    order = np.argsort(-scores, kind="stable")
    cand: dict[int, list[tuple[float, int]]] = {}
    for a, b, v in zip(pp, gg, iou):
        if v >= t:
            cand.setdefault(int(a), []).append((float(v), int(b)))
    taken, tp = set(), np.zeros(len(order))
    for k, p in enumerate(order):
        for v, g in sorted(cand.get(int(p), []), key=lambda x: (-x[0], x[1])):
            if g not in taken:
                taken.add(g)
                tp[k] = 1
                break
    ctp = np.cumsum(tp)
    precision = ctp / np.arange(1, len(order) + 1)
    recall = ctp / n_gt
    precision = np.maximum.accumulate(precision[::-1])[::-1] if len(precision) else precision
    points = []
    for r in np.linspace(0, 1, 101):
        idx = np.searchsorted(recall, r, side="left")
        points.append(precision[idx] if idx < len(precision) else 0.0)
    return float(np.mean(points))


# --------------------------------------------------------------------------- morphology
def _descriptors(mask: np.ndarray, offset, voxel_size: float) -> dict[str, float]:
    from skimage.measure import marching_cubes, mesh_surface_area

    coords = np.argwhere(mask)
    volume = float(len(coords)) * voxel_size ** 3
    centroid = (coords.mean(axis=0) + np.asarray(offset)) * voxel_size
    padded = np.pad(mask, 1).astype(np.float32)
    verts, faces, _, _ = marching_cubes(padded, level=0.5, spacing=(voxel_size,) * 3)
    surface = float(mesh_surface_area(verts, faces))
    d_eq = (6.0 * volume / np.pi) ** (1.0 / 3.0)
    sphericity = (np.pi ** (1 / 3) * (6 * volume) ** (2 / 3) / surface) if surface > 0 else float("nan")
    centred = (coords - coords.mean(axis=0)) * voxel_size
    if len(coords) > 1:
        eig = np.sort(np.clip(np.linalg.eigvalsh(np.cov(centred.T, bias=True)), 0, None))
        elongation = float(np.sqrt(eig[-1] / eig[0])) if eig[0] > 0 else float("inf")
    else:
        elongation = float("nan")
    return {"volume": volume, "d_eq": d_eq, "centroid_z": centroid[0], "centroid_y": centroid[1],
            "centroid_x": centroid[2], "surface_area": surface, "sphericity": sphericity,
            "elongation": elongation}


def morphology_table(pred: np.ndarray, gt: np.ndarray, pairs: list[tuple[int, int, float]],
                     voxel_size: float = 1.0) -> pd.DataFrame:
    """pairs: (prediction label, original GT label, IoU). Full objects, not the evaluation domain."""
    rows = []
    if not pairs:
        return pd.DataFrame(rows)
    pred_boxes = ndi.find_objects(pred.astype(np.int64, copy=False), max_label=max(p for p, _, _ in pairs))
    gt_boxes = ndi.find_objects(gt.astype(np.int64, copy=False), max_label=max(g for _, g, _ in pairs))
    for p, g, v in pairs:
        dp = _descriptors(pred[pred_boxes[p - 1]] == p, [s.start for s in pred_boxes[p - 1]], voxel_size)
        dg = _descriptors(gt[gt_boxes[g - 1]] == g, [s.start for s in gt_boxes[g - 1]], voxel_size)
        row = {"prediction_label": p, "gt_label": g, "iou": v}
        row |= {f"gt_{k}": x for k, x in dg.items()} | {f"pred_{k}": x for k, x in dp.items()}
        row["volume_rel_error"] = (dp["volume"] - dg["volume"]) / dg["volume"]
        row["d_eq_error"] = dp["d_eq"] - dg["d_eq"]
        row["centroid_error"] = float(np.linalg.norm([dp[f"centroid_{a}"] - dg[f"centroid_{a}"] for a in "zyx"]))
        row["surface_area_rel_error"] = (dp["surface_area"] - dg["surface_area"]) / dg["surface_area"]
        row["sphericity_error"] = dp["sphericity"] - dg["sphericity"]
        row["elongation_error"] = dp["elongation"] - dg["elongation"]
        rows.append(row)
    return pd.DataFrame(rows)


def _morph_summary(m: pd.DataFrame) -> dict[str, float]:
    if m.empty:
        return {"morph_n_matched": 0}
    out = {"morph_n_matched": len(m)}
    for col in ("volume_rel_error", "d_eq_error", "surface_area_rel_error", "sphericity_error", "elongation_error"):
        x = m[col].replace([np.inf, -np.inf], np.nan).dropna()
        out[f"{col}_mean"] = float(x.mean()) if len(x) else float("nan")
        out[f"{col}_median_abs"] = float(x.abs().median()) if len(x) else float("nan")
    out["centroid_error_mean"] = float(m["centroid_error"].mean())
    out["centroid_error_median"] = float(m["centroid_error"].median())
    return out


# --------------------------------------------------------------------------- driver
def evaluate(pred: np.ndarray, gt: np.ndarray, markers: np.ndarray | None = None, *,
             rule: str | None = None, ignore_overlap_threshold: float = 0.50, min_gt_voxels: int = 1000,
             exclude_face_touching: bool = True,
             thresholds=THRESHOLDS, scores: dict[int, float] | None = None, morph_iou: float | None = 0.50,
             voxel_size: float = 1.0) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """Returns (metrics row, per-prediction ignore table, morphology table of matched objects)."""
    pred = np.asarray(pred)
    if (pred < 0).any():
        raise ValueError("prediction contains negative labels")
    pred = pred.astype(np.int64, copy=False)
    protocol = reference_protocol(gt, min_gt_voxels, exclude_face_touching)
    table = classify_predictions(pred, protocol, markers, rule, ignore_overlap_threshold)
    eval_labels = np.sort(table.loc[~table["ignored"], "prediction_label"].to_numpy().astype(np.int64))
    # a prediction lying entirely inside excluded GT has no voxel in the domain: not evaluable (as e09b)
    domain_labels = np.unique(pred[~protocol["ignore_mask"]])
    eval_labels = eval_labels[np.isin(eval_labels, domain_labels)]
    pp, gg, iou, n_pred, n_gt = _iou_edges(pred, protocol, eval_labels)
    cat = table["support_category"].value_counts()
    out: dict[str, Any] = {
        "support": table["support"].iloc[0] if len(table) else ("seed_region" if markers is not None else "own_voxels"),
        "rule": table["rule"].iloc[0] if len(table) else rule,
        "ignore_overlap_threshold": ignore_overlap_threshold,
        "n_pred_total": int(len(table)), "n_pred_ignored": int(table["ignored"].sum()),
        "n_pred_background_support": int(cat.get("background", 0)),
        "n_pred_unseeded": int(cat.get("unseeded", 0)),
        "prediction_count_evaluable": n_pred, "valid_gt_count": n_gt,
        "n_gt_excluded": int(len(protocol["excluded_ids"])),
    }
    score_vec = None
    if scores is not None:
        score_vec = np.asarray([scores.get(int(p), np.nan) for p in eval_labels], dtype=float)
        if np.isnan(score_vec).any():
            raise ValueError("scores missing for some evaluated prediction labels")
    morph = pd.DataFrame()
    for t in thresholds:
        s = _suffix(t)
        rp, cg, miou = _max_card_matching(pp, gg, iou, n_pred, n_gt, t)
        tp = len(miou)
        fp, fn = n_pred - tp, n_gt - tp
        precision = tp / n_pred if n_pred else 0.0
        recall = tp / n_gt if n_gt else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        out |= {f"tp_iou{s}": tp, f"fp_iou{s}": fp, f"fn_iou{s}": fn, f"precision_iou{s}": precision,
                f"recall_iou{s}": recall, f"f1_iou{s}": f1,
                f"mean_matched_iou{s}": float(np.mean(np.sort(miou))) if tp else 0.0,
                f"ap_iou{s}": tp / (tp + fp + fn) if tp + fp + fn else 0.0,
                f"coco_ap_iou{s}": coco_ap(pp, gg, iou, score_vec, n_gt, t) if score_vec is not None else np.nan}
        if morph_iou is not None and abs(t - morph_iou) < 1e-9:
            pairs = [(int(eval_labels[a]), int(protocol["valid_ids"][b]), float(v)) for a, b, v in zip(rp, cg, miou)]
            morph = morphology_table(pred, protocol["gt"], pairs, voxel_size)
            out |= {"morph_iou": morph_iou, **_morph_summary(morph)}
    out["pred_to_gt_ratio"] = n_pred / n_gt if n_gt else float("nan")
    return out, table, morph


def main(argv=None) -> int:
    from .common import read3d

    ap = argparse.ArgumentParser(description="Label-agnostic seed/overlap-ignore evaluator (paper Sec. III-C)")
    ap.add_argument("--pred", type=Path, required=True, help="predicted instance labels (3D TIFF)")
    ap.add_argument("--gt", type=Path, required=True, help="reference instance labels (3D TIFF)")
    ap.add_argument("--markers", type=Path, help="markers that generated the prediction (optional)")
    ap.add_argument("--rule", choices=["dominant", "fraction"],
                    help="default: dominant with --markers, fraction without")
    ap.add_argument("--ignore-overlap-threshold", type=float, default=0.50)
    ap.add_argument("--min-gt-voxels", type=int, default=1000)
    ap.add_argument("--scores", type=Path, help="CSV prediction_label,score -> COCO-style AP")
    ap.add_argument("--morph-iou", type=float, default=0.50)
    ap.add_argument("--voxel-size", type=float, default=1.0, help="isotropic voxel size (output units)")
    ap.add_argument("--crop-id", default="")
    ap.add_argument("--method", default="")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    scores = None
    if args.scores:
        s = pd.read_csv(args.scores)
        scores = dict(zip(s["prediction_label"].astype(int), s["score"].astype(float)))
    markers = read3d(args.markers, "markers") if args.markers else None
    metrics, table, morph = evaluate(read3d(args.pred, "prediction"), read3d(args.gt, "GT"), markers,
                                     rule=args.rule, ignore_overlap_threshold=args.ignore_overlap_threshold,
                                     min_gt_voxels=args.min_gt_voxels, scores=scores, morph_iou=args.morph_iou,
                                     voxel_size=args.voxel_size)
    head = {"crop_id": args.crop_id, "method": args.method, "pred": str(args.pred), "gt": str(args.gt),
            "markers": str(args.markers) if args.markers else "", "voxel_size": args.voxel_size}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([head | metrics]).to_csv(args.out_dir / "metrics.csv", index=False)
    table.to_csv(args.out_dir / "prediction_ignore_table.csv", index=False)
    morph.to_csv(args.out_dir / "morphology_matched.csv", index=False)
    (args.out_dir / "metrics.json").write_text(json.dumps(head | metrics, indent=2, default=float), encoding="utf-8")
    print(json.dumps({k: metrics[k] for k in ("n_pred_ignored", "prediction_count_evaluable", "valid_gt_count",
                                              "f1_iou025", "f1_iou050", "f1_iou075")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
