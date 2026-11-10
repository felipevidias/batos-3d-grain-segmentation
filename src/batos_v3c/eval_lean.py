"""Memory-lean, exact re-statement of batos_v2.evaluate_v2.evaluate (metrics row only; morph_iou=None, scores=None).

Why: the frozen evaluator casts whole volumes to int64 and builds several full-size temporaries; on the 400^3 / 420^3
canonical-field-of-view windows (Phase 5i) one call peaks at ~4-5 GB, above the 4 GB working rule. Every quantity it uses
is a voxel count -- per prediction label, per (prediction, GT) pair on the support, per (prediction, valid GT) pair on
the domain ~ignore_mask -- so the counts are accumulated over z-slabs here and the frozen downstream code (ignore rule,
_max_card_matching, metric formulas) is then applied unchanged. The frozen module is imported, never modified.
Equality with the frozen evaluator is checked in scripts/check_eval_lean.py (all metric keys, exact).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from batos_v2.evaluate_v2 import THRESHOLDS, _max_card_matching, _suffix

SLAB = 24


def _slabs(n: int, step: int = SLAB):
    for z in range(0, n, step):
        yield slice(z, min(n, z + step))


def _merge(parts: list[pd.Series]) -> pd.Series:
    if not parts:
        return pd.Series(dtype=np.int64)
    return pd.concat(parts).groupby(level=list(range(parts[0].index.nlevels))).sum().astype(np.int64)


def protocol(gt: np.ndarray, min_gt_voxels: int = 1000, exclude_face_touching: bool = True) -> dict:
    """reference_protocol without full-size arrays: areas, valid / excluded ids, lookup, remap."""
    if gt.ndim != 3:
        raise ValueError(f"GT must be 3D, got {gt.shape}")
    if np.issubdtype(gt.dtype, np.signedinteger) and (gt < 0).any():
        raise ValueError("GT contains negative labels")
    G = int(gt.max()) + 1
    areas = np.zeros(G, np.int64)
    for s in _slabs(gt.shape[0]):
        areas += np.bincount(gt[s].ravel(), minlength=G)
    border = set()
    if exclude_face_touching:
        faces = np.concatenate([gt[0].ravel(), gt[-1].ravel(), gt[:, 0].ravel(), gt[:, -1].ravel(),
                                gt[:, :, 0].ravel(), gt[:, :, -1].ravel()]).astype(np.int64)
        border = {int(v) for v in np.unique(faces) if v > 0}
    valid_ids = np.asarray([lab for lab in range(1, len(areas))
                            if areas[lab] >= min_gt_voxels and lab not in border], dtype=np.int64)
    present = np.flatnonzero(areas)
    vset = set(valid_ids.tolist())
    excluded_ids = np.asarray([lab for lab in present if lab > 0 and lab not in vset], dtype=np.int64)
    valid_lookup = np.zeros(len(areas), dtype=bool)
    valid_lookup[valid_ids] = True
    remap = np.zeros(len(areas), dtype=np.int64)
    remap[valid_ids] = np.arange(1, len(valid_ids) + 1)
    return {"areas": areas, "valid_ids": valid_ids, "excluded_ids": excluded_ids, "valid_lookup": valid_lookup,
            "remap": remap}


def evaluate(pred: np.ndarray, gt: np.ndarray, markers: np.ndarray | None = None, *, rule: str | None = None,
             ignore_overlap_threshold: float = 0.50, min_gt_voxels: int = 1000, exclude_face_touching: bool = True,
             thresholds=THRESHOLDS, proto: dict | None = None) -> dict[str, Any]:
    """Same metrics row as batos_v2.evaluate_v2.evaluate(..., morph_iou=None, scores=None)[0]."""
    if rule is None:
        rule = "dominant" if markers is not None else "fraction"
    if rule not in ("dominant", "fraction"):
        raise ValueError(f"unknown rule {rule!r}")
    if pred.shape != gt.shape or (markers is not None and markers.shape != gt.shape):
        raise ValueError(f"shape mismatch pred={pred.shape} gt={gt.shape} "
                         f"markers={None if markers is None else markers.shape}")
    if np.issubdtype(pred.dtype, np.signedinteger) and (pred < 0).any():
        raise ValueError("prediction contains negative labels")
    proto = proto or protocol(gt, min_gt_voxels, exclude_face_touching)
    lookup, remap = proto["valid_lookup"], proto["remap"]
    n_gt = len(proto["valid_ids"])
    vox, sup, pos_, exc, pairs, joint, dom_lab = [], [], [], [], [], [], []
    for s in _slabs(gt.shape[0]):
        p = pred[s].astype(np.int64).ravel()
        g = gt[s].astype(np.int64).ravel()
        vg = lookup[g]
        ign = (g > 0) & ~vg
        pp = p[p > 0]
        vox.append(pd.Series(pp).value_counts())
        support = (p > 0) if markers is None else ((p > 0) & (markers[s].ravel() > 0))
        sp, sg, sx = p[support], g[support], ign[support]
        sup.append(pd.Series(sp).value_counts())
        pos_.append(pd.Series(sp[sg > 0]).value_counts())
        exc.append(pd.Series(sp[sx]).value_counts())
        k = sg > 0
        if k.any():
            pairs.append(pd.DataFrame({"p": sp[k], "g": sg[k]}).value_counts())
        dom = ~ign
        pv, gv = p[dom], np.where(vg, remap[g], 0)[dom]
        dom_lab.append(np.unique(pv))
        joint.append(pd.DataFrame({"p": pv, "g": gv}).value_counts())
    # ---- classify_predictions (frozen logic on aggregated counts)
    vox_s = _merge(vox).sort_index()
    table = pd.DataFrame({"prediction_label": vox_s.index.to_numpy(np.int64), "prediction_voxels": vox_s.to_numpy(np.int64)})
    table["support_voxels"] = table["prediction_label"].map(_merge(sup)).fillna(0).astype(np.int64)
    table["support_gt_positive_voxels"] = table["prediction_label"].map(_merge(pos_)).fillna(0).astype(np.int64)
    table["support_excluded_voxels"] = table["prediction_label"].map(_merge(exc)).fillna(0).astype(np.int64)
    table["support_excluded_fraction"] = (table["support_excluded_voxels"]
                                          / table["support_voxels"].where(table["support_voxels"] > 0))
    if pairs:
        pr = _merge(pairs).rename("n").reset_index()
        pr = pr.sort_values(["p", "n", "g"], ascending=[True, False, True]).drop_duplicates("p")
        dominant = pr.set_index("p")["g"]
    else:
        dominant = pd.Series(dtype=np.int64)
    table["dominant_gt_label"] = table["prediction_label"].map(dominant).fillna(0).astype(np.int64)
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
        ignored = (table["support_excluded_fraction"].fillna(0.0) >= ignore_overlap_threshold).to_numpy()
    table["ignored"] = ignored
    # ---- evaluate / _iou_edges (frozen logic on aggregated counts)
    eval_labels = np.sort(table.loc[~table["ignored"], "prediction_label"].to_numpy().astype(np.int64))
    domain_labels = np.unique(np.concatenate(dom_lab)) if dom_lab else np.zeros(0, np.int64)
    eval_labels = eval_labels[np.isin(eval_labels, domain_labels)]
    j = _merge(joint)
    jp = j.index.get_level_values("p").to_numpy(np.int64) if len(j) else np.zeros(0, np.int64)
    jg = j.index.get_level_values("g").to_numpy(np.int64) if len(j) else np.zeros(0, np.int64)
    jn = j.to_numpy(np.int64)
    n_pred = len(eval_labels)
    index = np.zeros(int(pred.max()) + 1 if pred.size else 1, dtype=np.int64)
    index[eval_labels] = np.arange(1, n_pred + 1)
    pc = index[jp]
    pred_areas = np.bincount(pc, weights=jn, minlength=n_pred + 1).astype(np.int64)
    gt_areas = np.bincount(jg, weights=jn, minlength=n_gt + 1).astype(np.int64)
    both = (pc > 0) & (jg > 0)
    code = pc[both] * (n_gt + 1) + jg[both]
    order = np.argsort(code, kind="stable")
    code, inter = code[order], jn[both][order]
    pp, gg = code // (n_gt + 1), code % (n_gt + 1)
    iou = inter / np.maximum(pred_areas[pp] + gt_areas[gg] - inter, 1)
    pp, gg = pp - 1, gg - 1
    cat = table["support_category"].value_counts()
    out: dict[str, Any] = {
        "support": "seed_region" if markers is not None else "own_voxels", "rule": rule,
        "ignore_overlap_threshold": ignore_overlap_threshold,
        "n_pred_total": int(len(table)), "n_pred_ignored": int(table["ignored"].sum()),
        "n_pred_background_support": int(cat.get("background", 0)),
        "n_pred_unseeded": int(cat.get("unseeded", 0)),
        "prediction_count_evaluable": n_pred, "valid_gt_count": n_gt,
        "n_gt_excluded": int(len(proto["excluded_ids"])),
    }
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
    out["pred_to_gt_ratio"] = n_pred / n_gt if n_gt else float("nan")
    return out


def matched(pred: np.ndarray, gt: np.ndarray, valid) -> dict[int, tuple[int, float, float]]:
    """valid GT id -> (pred id, IoU, relative volume error) for IoU > 0.5 (slab-wise; same as run_v3c_dev.matched)."""
    P, G = int(pred.max()) + 1, int(gt.max()) + 1
    sp, sg, parts = np.zeros(P, np.int64), np.zeros(G, np.int64), []
    for s in _slabs(gt.shape[0]):
        p, g = pred[s].astype(np.int64).ravel(), gt[s].astype(np.int64).ravel()
        sp += np.bincount(p, minlength=P)
        sg += np.bincount(g, minlength=G)
        b = (p > 0) & (g > 0)
        parts.append(pd.DataFrame({"p": p[b], "g": g[b]}).value_counts())
    j = _merge(parts)
    if not len(j):
        return {}
    pu = j.index.get_level_values("p").to_numpy(np.int64)
    gu = j.index.get_level_values("g").to_numpy(np.int64)
    inter = j.to_numpy(np.int64)
    iou = inter / (sp[pu] + sg[gu] - inter)
    vs = set(np.asarray(valid).tolist())
    return {int(b): (int(a), float(i), float(sp[a] / sg[b] - 1)) for a, b, i in zip(pu, gu, iou)
            if i > 0.5 and int(b) in vs}


def seed_stats(markers: np.ndarray, gt: np.ndarray, proto: dict) -> dict:
    """run_v32_dev.seed_stats with the GT areas / valid ids taken from protocol() (no full-size int64 copy)."""
    pos = np.flatnonzero(markers.ravel())
    mk, g = markers.ravel()[pos].astype(np.int64), gt.ravel()[pos].astype(np.int64)
    df = pd.DataFrame({"m": mk, "g": g}).value_counts().reset_index(name="n").sort_values(["m", "n"],
                                                                                          ascending=[True, False])
    owner = df.drop_duplicates("m")
    per = owner[owner.g > 0].g.value_counts()
    valid = proto["valid_ids"]
    cnt = per.reindex(valid, fill_value=0).to_numpy()
    n_all = int(len(np.flatnonzero(proto["areas"][1:])))
    return {"n_markers": int(len(owner)), "markers_background": int((owner.g == 0).sum()), "gt_all": n_all,
            "gt_valid": int(len(valid)), "markers_per_grain_all": len(owner) / n_all if n_all else None,
            "markers_per_valid_grain": float(cnt.mean()) if len(cnt) else None,
            "valid_seeded": int((cnt >= 1).sum()), "valid_unseeded": int((cnt == 0).sum()),
            "valid_split": int((cnt >= 2).sum()),
            "frac_seeded": float((cnt >= 1).mean()) if len(cnt) else None,
            "frac_split": float((cnt >= 2).mean()) if len(cnt) else None}
