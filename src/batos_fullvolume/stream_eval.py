"""Streaming (z-slab) evaluator equivalent to batos_v2.evaluate_v2.evaluate (rule "dominant", markers).

Only the accumulation is re-implemented: per-label areas and (prediction, valid GT) intersections are summed over
z-slabs, the seed-ignore decision is taken from the seed voxels, and matching / metrics call evaluate_v2's own
_max_card_matching and formulas. The reference protocol is given as an explicit set of valid GT ids instead of
evaluate_v2's "touches no crop face" rule (full volume: ROI list + >= 1000 voxels); every other GT label present is
excluded (ignore mask), exactly as in evaluate_v2.reference_protocol. scripts/verify_stream_eval.py checks that
this module reproduces evaluate_v2.evaluate on dev crops (face rule, small slabs) key by key.

Extra outputs (not in evaluate_v2, defined here a priori):
  overseg_ratio = n_pred_eval / n_valid_gt (repo convention);
  split GT  = valid GT g with >= 2 evaluated predictions p such that |p n g| >= 0.5 |p| (domain voxels);
  merge prediction = evaluated p with >= 2 valid GT g such that |p n g| >= 0.5 |g|.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.ndimage as ndi

from batos_v2.evaluate_v2 import THRESHOLDS, _max_card_matching, _suffix


def gt_label_stats(gt, slab: int = 16, box=None) -> pd.DataFrame:
    """Per GT label (over `box`, default whole array): voxels, centroid, bbox (global coordinates)."""
    box = box or tuple((0, s) for s in gt.shape)
    (z0, z1), (y0, y1), (x0, x1) = box
    n = int(np.iinfo(gt.dtype).max) + 1 if gt.dtype.kind == "u" and gt.dtype.itemsize <= 2 else None
    acc = None
    for zs in range(z0, z1, slab):
        ze = min(z1, zs + slab)
        g = np.asarray(gt[zs:ze, y0:y1, x0:x1])
        if n is None:
            n = int(max(g.max(), 0)) + 1  # only reached for wide dtypes; grows below if needed
        idx = np.flatnonzero(g)
        if acc is None:
            acc = {k: np.zeros(n, np.float64) for k in ("v", "sz", "sy", "sx")}
            acc |= {"lo": np.full((n, 3), np.iinfo(np.int64).max), "hi": np.full((n, 3), -1)}
        lab = g.ravel()[idx].astype(np.int64)
        zz, yy, xx = np.unravel_index(idx, g.shape)
        zz = zz + zs
        yy = yy + y0
        xx = xx + x0
        acc["v"] += np.bincount(lab, minlength=n)
        for k, c in (("sz", zz), ("sy", yy), ("sx", xx)):
            acc[k] += np.bincount(lab, weights=c, minlength=n)
        for i, sl in enumerate(ndi.find_objects(g)):
            if sl is None:
                continue
            lo = (sl[0].start + zs, sl[1].start + y0, sl[2].start + x0)
            hi = (sl[0].stop - 1 + zs, sl[1].stop - 1 + y0, sl[2].stop - 1 + x0)
            acc["lo"][i + 1] = np.minimum(acc["lo"][i + 1], lo)
            acc["hi"][i + 1] = np.maximum(acc["hi"][i + 1], hi)
    ids = np.flatnonzero(acc["v"])
    ids = ids[ids > 0]
    v = acc["v"][ids]
    return pd.DataFrame({"gt_label": ids, "voxels": v.astype(np.int64), "cz": acc["sz"][ids] / v,
                         "cy": acc["sy"][ids] / v, "cx": acc["sx"][ids] / v,
                         "z0": acc["lo"][ids, 0], "y0": acc["lo"][ids, 1], "x0": acc["lo"][ids, 2],
                         "z1": acc["hi"][ids, 0] + 1, "y1": acc["hi"][ids, 1] + 1, "x1": acc["hi"][ids, 2] + 1})


def accumulate(pred, gt, valid_ids: np.ndarray, pred_offset=(0, 0, 0), slab: int = 16) -> dict:
    """One streaming pass. pred: array-like (memmap) over a box starting at pred_offset in GT coordinates;
    gt: array-like over the full volume (memmap). Returns raw counts on the pred box."""
    o = np.asarray(pred_offset)
    zb, yb, xb = pred.shape
    gmax = int(np.iinfo(gt.dtype).max) if gt.dtype.kind == "u" else int(valid_ids.max(initial=0))
    valid_lookup = np.zeros(gmax + 1, bool)
    valid_lookup[valid_ids] = True
    remap = np.zeros(gmax + 1, np.int64)
    remap[valid_ids] = np.arange(1, len(valid_ids) + 1)
    n_gt = len(valid_ids)
    pmax = 0
    pred_total, pred_dom, gt_area = np.zeros(1, np.int64), np.zeros(1, np.int64), np.zeros(n_gt + 1, np.int64)
    inter_parts = []
    gt_present = np.zeros(gmax + 1, bool)
    for zs in range(0, zb, slab):
        ze = min(zb, zs + slab)
        p = np.asarray(pred[zs:ze]).astype(np.int64, copy=False)
        g = np.asarray(gt[o[0] + zs:o[0] + ze, o[1]:o[1] + yb, o[2]:o[2] + xb]).astype(np.int64, copy=False)
        gt_present[np.unique(g)] = True
        pm = int(p.max(initial=0))
        if pm + 1 > len(pred_total):
            pred_total = np.pad(pred_total, (0, pm + 1 - len(pred_total)))
            pred_dom = np.pad(pred_dom, (0, pm + 1 - len(pred_dom)))
        pmax = max(pmax, pm)
        pred_total += np.bincount(p.ravel(), minlength=len(pred_total))
        dom = ~((g > 0) & ~valid_lookup[g])
        pv = p[dom]
        gv = remap[g[dom]]
        pred_dom += np.bincount(pv, minlength=len(pred_dom))
        gt_area += np.bincount(gv, minlength=n_gt + 1)
        both = (pv > 0) & (gv > 0)
        if both.any():
            code, cnt = np.unique(pv[both] * (n_gt + 1) + gv[both], return_counts=True)
            inter_parts.append(pd.DataFrame({"code": code, "n": cnt}))
    inter = (pd.concat(inter_parts).groupby("code", sort=True)["n"].sum() if inter_parts
             else pd.Series(dtype=np.int64))
    codes = inter.index.to_numpy(np.int64)
    return {"pred_total": pred_total, "pred_dom": pred_dom, "gt_area": gt_area, "n_gt": n_gt,
            "valid_ids": np.asarray(valid_ids, np.int64), "valid_lookup": valid_lookup,
            "gt_present": np.flatnonzero(gt_present[1:]) + 1,
            "inter_p": codes // (n_gt + 1), "inter_g": codes % (n_gt + 1), "inter_n": inter.to_numpy(np.int64)}


def seed_values(pred, gt, seeds_zyx: np.ndarray, pred_offset=(0, 0, 0), chunk: int = 2_000_000):
    """(pred value, GT value) at every seed voxel (global zyx) that lies in the pred box."""
    o = np.asarray(pred_offset)
    loc = seeds_zyx - o
    ok = np.all((loc >= 0) & (loc < np.asarray(pred.shape)), axis=1)
    loc, glob = loc[ok], seeds_zyx[ok]
    order = np.lexsort((loc[:, 2], loc[:, 1], loc[:, 0]))
    loc, glob = loc[order], glob[order]
    pv = np.empty(len(loc), np.int64)
    gv = np.empty(len(loc), np.int64)
    for s in range(0, len(loc), chunk):  # slab-wise reads keep memmap access sequential
        zl = loc[s:s + chunk]
        zg = glob[s:s + chunk]
        for z in np.unique(zl[:, 0]):
            m = zl[:, 0] == z
            row_p = np.asarray(pred[int(z)])
            row_g = np.asarray(gt[int(zg[m][0, 0])])
            pv[s:s + chunk][m] = row_p[zl[m, 1], zl[m, 2]]
            gv[s:s + chunk][m] = row_g[zg[m, 1], zg[m, 2]]
    return pv, gv


def finalize(acc: dict, seed_pred: np.ndarray, seed_gt: np.ndarray, remove: np.ndarray | None = None,
             thresholds=THRESHOLDS) -> tuple[dict, pd.DataFrame, dict]:
    """evaluate_v2.evaluate (rule dominant, morph off) from accumulated counts. `remove` = prediction labels
    deleted by a size filter (they vanish from the segmentation, as size_filter() in run_test_v3)."""
    pred_total = acc["pred_total"]
    pred_ids = np.flatnonzero(pred_total)
    pred_ids = pred_ids[pred_ids > 0]
    if remove is not None and len(remove):
        pred_ids = pred_ids[~np.isin(pred_ids, remove)]
    lookup = acc["valid_lookup"]
    # --- classify_predictions (markers, rule dominant)
    m = (seed_pred > 0) & np.isin(seed_pred, pred_ids)
    sp, sg = seed_pred[m], seed_gt[m]
    s_ex = (sg > 0) & ~lookup[sg]
    table = pd.DataFrame({"prediction_label": pred_ids, "prediction_voxels": pred_total[pred_ids]})
    n_sup = pd.Series(sp).value_counts()
    n_pos = pd.Series(sp[sg > 0]).value_counts()
    n_exc = pd.Series(sp[s_ex]).value_counts()
    table["support_voxels"] = table["prediction_label"].map(n_sup).fillna(0).astype(np.int64)
    table["support_gt_positive_voxels"] = table["prediction_label"].map(n_pos).fillna(0).astype(np.int64)
    table["support_excluded_voxels"] = table["prediction_label"].map(n_exc).fillna(0).astype(np.int64)
    pos = sg > 0
    if pos.any():
        pairs = pd.DataFrame({"p": sp[pos], "g": sg[pos]}).value_counts().rename("n").reset_index()
        pairs = pairs.sort_values(["p", "n", "g"], ascending=[True, False, True]).drop_duplicates("p")
        dominant = pairs.set_index("p")["g"]
    else:
        dominant = pd.Series(dtype=np.int64)
    table["dominant_gt_label"] = table["prediction_label"].map(dominant).fillna(0).astype(np.int64)
    dom = table["dominant_gt_label"].to_numpy()
    dom_valid = np.zeros(len(dom), bool)
    inside = (dom > 0) & (dom < len(lookup))
    dom_valid[inside] = lookup[dom[inside]]
    table["support_category"] = np.where(table["support_voxels"] == 0, "unseeded",
                                         np.where(dom == 0, "background",
                                                  np.where(dom_valid, "valid_gt", "ignored_gt")))
    table["ignored"] = table["support_category"] == "ignored_gt"
    # --- evaluated predictions: not ignored and >= 1 voxel in the domain
    eval_labels = np.sort(table.loc[~table["ignored"], "prediction_label"].to_numpy().astype(np.int64))
    pd_ = acc["pred_dom"]
    eval_labels = eval_labels[(eval_labels < len(pd_)) & (pd_[np.minimum(eval_labels, len(pd_) - 1)] > 0)]
    n_pred, n_gt = len(eval_labels), acc["n_gt"]
    index = np.zeros(len(pred_total), np.int64)
    index[eval_labels] = np.arange(1, n_pred + 1)
    pc = index[acc["inter_p"]]
    keep = pc > 0
    pp, gg, inter = pc[keep], acc["inter_g"][keep], acc["inter_n"][keep]
    pa = pd_[eval_labels]
    ga = acc["gt_area"][1:]
    iou = inter / np.maximum(pa[pp - 1] + ga[gg - 1] - inter, 1)
    pp, gg = pp - 1, gg - 1
    cat = table["support_category"].value_counts()
    out = {"support": "seed_region", "rule": "dominant", "ignore_overlap_threshold": 0.5,
           "n_pred_total": int(len(table)), "n_pred_ignored": int(table["ignored"].sum()),
           "n_pred_background_support": int(cat.get("background", 0)),
           "n_pred_unseeded": int(cat.get("unseeded", 0)),
           "prediction_count_evaluable": n_pred, "valid_gt_count": n_gt,
           "n_gt_excluded": int(len(np.setdiff1d(acc["gt_present"], acc["valid_ids"])))}
    matches = {}
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
                f"ap_iou{s}": tp / (tp + fp + fn) if tp + fp + fn else 0.0, f"coco_ap_iou{s}": np.nan}
        matches[t] = (eval_labels[rp], acc["valid_ids"][cg], miou)
    out["pred_to_gt_ratio"] = n_pred / n_gt if n_gt else float("nan")
    # --- over / under-segmentation (extra, defined in the module docstring)
    if len(pp):
        frac_p = inter / pa[pp]
        frac_g = inter / ga[gg]
        split = pd.Series(gg[frac_p >= 0.5]).value_counts()
        merge = pd.Series(pp[frac_g >= 0.5]).value_counts()
        n_split, n_merge = int((split >= 2).sum()), int((merge >= 2).sum())
        n_gt_in_merge = int(merge[merge >= 2].sum())
    else:
        n_split = n_merge = n_gt_in_merge = 0
    out |= {"overseg_ratio": n_pred / n_gt if n_gt else float("nan"), "n_gt_split": n_split,
            "frac_gt_split": n_split / n_gt if n_gt else float("nan"), "n_pred_merge": n_merge,
            "n_gt_in_merge": n_gt_in_merge, "frac_gt_in_merge": n_gt_in_merge / n_gt if n_gt else float("nan")}
    detail = {"eval_labels": eval_labels, "matches": matches, "pp": pp, "gg": gg, "iou": iou}
    return out, table, detail
