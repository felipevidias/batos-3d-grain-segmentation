#!/usr/bin/env python3
"""Phase 2: shape fidelity, boundary accuracy and grain contacts -- cpsam_v2 vs BA-TOS v3 vs GT (dev + test).

Per crop (64), per method (cellpose:cpsam_v2, v3 = BA-TOS_sel), size filter off:
1. Matching: evaluate_v2 no-marker mode (fraction 0.50), pairs at IoU >= 0.50 (unique at t >= 0.5).
   Morphology per matched pair from evaluate_v2.morphology_table (full objects, voxel units): volume, d_eq,
   surface area (marching cubes), sphericity, elongation, centroid, and their errors (pred - GT).
2. Boundary metrics per matched pair: surface voxels = mask minus its 6-conn erosion; distances by EDT;
   MSD = mean of both directed surface distances pooled; HD95 = max of the two directed 95th percentiles.
3. Alternative scoring modes (summary only): strict_noignore (nothing ignored, standard GT exclusion),
   allgt_faces (nothing ignored, face-touching GT >= 1000 voxels counted), allgt_all (nothing ignored, every GT
   label counted, min size 1).
4. Contacts: 6-connectivity adjacency pairs (two different labels on face-adjacent voxels). Reference = pairs of
   two VALID GT grains (>= 1000 voxels, not touching a face). recall = reference contacts whose two grains are both
   matched and whose matched predictions are adjacent; precision = adjacent prediction pairs with both members
   matched to valid GT whose GT pair is in contact. Also compared with the authors' contact list
   (PhDUpload/<SPEC>/Contacts/ContactPoint/<scan>.txt: grain A, grain B, Z Y X in full-volume label ids; crop GT keeps
   the full-volume ids), restricted to pairs of valid GT grains of the crop.
Outputs: results/cellpose/shape/{pairs.csv.gz, crop_summary.csv, contacts.csv, modes.csv, tables.md}.
Batos env: python -I scripts\\shape_contacts_cellpose_v3.py --workers 8
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

CP = Path(r"C:\BA_TOS_DATA\runs\cellpose\cpsam_v2")
V3_DEV = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test\runs\v3\x1")
V3_TEST = Path(r"C:\BA_TOS_DATA\runs\test_protocol\v3")
V3_SEG = Path("e23c") / "watershed" / "segmentations" / "batos_selective_distance_c1.tif"
CONTACTS = Path(r"C:\Users\USER\Downloads\PhDUpload\PhDUpload")
OUT = REPO / "results" / "cellpose" / "shape"
SIX = ndi.generate_binary_structure(3, 1)
MODES = {"no_markers_fraction050": dict(rule="fraction", ignore_overlap_threshold=0.50),
         "strict_noignore": dict(rule="fraction", ignore_overlap_threshold=1.01),
         "allgt_faces": dict(rule="fraction", ignore_overlap_threshold=1.01, exclude_face_touching=False),
         "allgt_all": dict(rule="fraction", ignore_overlap_threshold=1.01, exclude_face_touching=False,
                           min_gt_voxels=1)}


def crops() -> pd.DataFrame:
    dev = pd.read_csv(REPO / "data" / "metadata" / "crops_manifest.csv").assign(split="dev", specimen="EFR",
                                                                                test_class="dev")
    test = pd.read_csv(REPO / "data" / "metadata" / "crops_manifest_test.csv").assign(split="test")
    cols = ["crop_id", "split", "specimen", "scan", "test_class", "x", "y", "z", "size", "labels_path"]
    return pd.concat([dev[cols], test[cols]], ignore_index=True)


def adjacency(lab: np.ndarray) -> set[tuple[int, int]]:
    pairs = []
    for ax in range(3):
        a = np.moveaxis(lab, ax, 0)
        u, v = a[:-1], a[1:]
        m = (u > 0) & (v > 0) & (u != v)
        if m.any():
            p = np.stack([np.minimum(u[m], v[m]), np.maximum(u[m], v[m])], 1)
            pairs.append(np.unique(p, axis=0))
    if not pairs:
        return set()
    return set(map(tuple, np.unique(np.concatenate(pairs), axis=0).tolist()))


def valid_labels(gt: np.ndarray, min_vox: int = 1000) -> set[int]:
    ids, n = np.unique(gt[gt > 0], return_counts=True)
    face = set()
    for sl in (gt[0], gt[-1], gt[:, 0], gt[:, -1], gt[:, :, 0], gt[:, :, -1]):
        face |= set(np.unique(sl).tolist())
    return {int(i) for i, k in zip(ids, n) if k >= min_vox and int(i) not in face}


def boundary(pm: np.ndarray, gm: np.ndarray) -> tuple[float, float]:
    sp = pm & ~ndi.binary_erosion(pm, SIX)
    sg = gm & ~ndi.binary_erosion(gm, SIX)
    dg = ndi.distance_transform_edt(~sg)
    dp = ndi.distance_transform_edt(~sp)
    a, b = dg[sp], dp[sg]
    return float(np.concatenate([a, b]).mean()), float(max(np.percentile(a, 95), np.percentile(b, 95)))


@lru_cache(maxsize=None)
def author_contacts(scan: str) -> pd.DataFrame | None:
    p = CONTACTS / scan[:3] / "Contacts" / "ContactPoint" / f"{scan}.txt"
    if not p.exists():
        return None
    d = pd.read_csv(p, sep=r"\s+", comment="#", header=None, names=["a", "b", "z", "y", "x"])
    return d


def one_crop(c: dict) -> dict:
    from batos_v2.evaluate_v2 import evaluate
    gt = tifffile.imread(c["labels_path"]).astype(np.int64)
    valid = valid_labels(gt)
    gt_adj = adjacency(gt)
    ref = {p for p in gt_adj if p[0] in valid and p[1] in valid}
    preds = {"cellpose:cpsam_v2": CP / c["crop_id"] / "labels.tif",
             "v3": (V3_DEV if c["split"] == "dev" else V3_TEST) / c["crop_id"] / V3_SEG}
    base = {k: c[k] for k in ("crop_id", "split", "specimen", "scan", "test_class")}
    pair_rows, mode_rows, contact_rows = [], [], []

    # authors' contacts restricted to pairs of valid GT grains of this crop
    au = author_contacts(c["scan"])
    au_set = None
    if au is not None:
        sub = au[au.a.isin(valid) & au.b.isin(valid)]
        au_set = {(min(a, b), max(a, b)) for a, b in zip(sub.a, sub.b)}
        inside = ((sub.z >= c["z"]) & (sub.z < c["z"] + c["size"]) & (sub.y >= c["y"]) & (sub.y < c["y"] + c["size"])
                  & (sub.x >= c["x"]) & (sub.x < c["x"] + c["size"]))
        tp = len(ref & au_set)
        contact_rows.append(base | {"method": "GT (6-adj)", "reference": "authors", "n_ref": len(au_set),
                                    "n_pred_pairs": len(ref), "tp": tp,
                                    "recall": tp / len(au_set) if au_set else np.nan,
                                    "precision": tp / len(ref) if ref else np.nan,
                                    "authors_points_inside_crop": float(inside.mean()) if len(sub) else np.nan})

    for method, path in preds.items():
        pred = tifffile.imread(path).astype(np.int64)
        for mode, kw in MODES.items():
            s, _, morph = evaluate(pred, gt, None, morph_iou=0.50 if mode == "no_markers_fraction050" else None, **kw)
            mode_rows.append(base | {"method": method, "mode": mode,
                                     **{k: s[k] for k in s if any(k.startswith(p) for p in
                                        ("tp_", "fp_", "fn_", "f1_", "precision_", "recall_", "n_pred", "valid_gt",
                                         "n_gt", "prediction_count", "mean_matched"))}})
            if mode != "no_markers_fraction050":
                continue
            # boundary metrics per matched pair
            pb = ndi.find_objects(pred)
            gb = ndi.find_objects(gt)
            msd, hd = [], []
            for p, g in zip(morph.prediction_label, morph.gt_label):
                sl = tuple(slice(max(0, min(a.start, b.start) - 2), max(a.stop, b.stop) + 2)
                           for a, b in zip(pb[p - 1], gb[g - 1]))
                m1, m2 = boundary(pred[sl] == p, gt[sl] == g)
                msd.append(m1)
                hd.append(m2)
            morph = morph.assign(msd=msd, hd95=hd, **base, method=method)
            pair_rows.append(morph)
            # contacts
            p2g = dict(zip(morph.prediction_label.astype(int), morph.gt_label.astype(int)))
            g2p = {g: p for p, g in p2g.items()}
            padj = adjacency(pred)
            for rname, rset in (("GT 6-adj", ref), ("authors", au_set)):
                if rset is None:
                    continue
                rec = sum(1 for a, b in rset if a in g2p and b in g2p
                          and (min(g2p[a], g2p[b]), max(g2p[a], g2p[b])) in padj)
                cand = [(p2g[a], p2g[b]) for a, b in padj if a in p2g and b in p2g]
                prec = sum(1 for a, b in cand if (min(a, b), max(a, b)) in rset)
                contact_rows.append(base | {"method": method, "reference": rname, "n_ref": len(rset),
                                            "n_pred_pairs": len(cand), "tp": rec, "tp_precision": prec,
                                            "recall": rec / len(rset) if rset else np.nan,
                                            "precision": prec / len(cand) if cand else np.nan,
                                            "n_pred_adjacent_total": len(padj)})
    return {"pairs": pd.concat(pair_rows, ignore_index=True) if pair_rows else pd.DataFrame(),
            "modes": pd.DataFrame(mode_rows), "contacts": pd.DataFrame(contact_rows)}


def scopes(df):
    yield "dev (12 EFR)", df[df.split == "dev"]
    t = df[df.split == "test"]
    yield "test: all 52", t
    yield "test: EFR", t[t.specimen == "EFR"]
    yield "test: external", t[t.specimen != "EFR"]
    for s in sorted(t.specimen.unique()):
        yield f"test: {s}", t[t.specimen == s]


def tables(pairs: pd.DataFrame, modes: pd.DataFrame, contacts: pd.DataFrame) -> str:
    L = ["# Shape, boundary and contacts: cpsam_v2 vs BA-TOS v3 vs GT (size filter off)", ""]
    # paired comparison on GT grains matched by both methods
    pairs = pairs.assign(volume_err_pct=100 * pairs.volume_rel_error, surface_err_pct=100 * pairs.surface_area_rel_error)
    key = ["crop_id", "gt_label"]
    both = pairs.groupby(key).method.nunique()
    common = pairs.set_index(key).loc[both[both == 2].index].reset_index()
    cols = [("volume_err_pct", "volume err %"), ("d_eq_error", "d_eq err (vox)"), ("surface_err_pct", "surface err %"),
            ("sphericity_error", "sphericity err"), ("elongation_error", "elongation err"),
            ("centroid_error", "centroid err (vox)"), ("msd", "MSD (vox)"), ("hd95", "HD95 (vox)")]
    for title, data in (("all matched pairs of each method", pairs), ("GT grains matched by BOTH methods (paired)", common)):
        L += [f"## Morphology and boundary, {title}: mean signed error [median |error|]", "",
              "| scope | method | n pairs | " + " | ".join(n for _, n in cols) + " |", "|---|---|---|" + "---|" * len(cols)]
        for scope, d in scopes(data):
            for method, g in d.groupby("method"):
                cells = []
                for col, _ in cols:
                    x = g[col].replace([np.inf, -np.inf], np.nan).dropna()
                    cells.append(f"{x.mean():+.3f} [{x.abs().median():.3f}]" if col not in ("msd", "hd95", "centroid_error")
                                 else f"{x.mean():.3f} [{x.median():.3f}]")
                L.append(f"| {scope} | {method} | {len(g)} | " + " | ".join(cells) + " |")
        L.append("")
    L += ["## Mean descriptor values on the paired set (GT vs prediction)", "",
          "| scope | method | GT sphericity | pred sphericity | GT elongation | pred elongation | GT d_eq | pred d_eq |",
          "|---|---|---|---|---|---|---|---|"]
    for scope, d in scopes(common):
        for method, g in d.groupby("method"):
            L.append(f"| {scope} | {method} | {g.gt_sphericity.mean():.3f} | {g.pred_sphericity.mean():.3f} | "
                     f"{g.gt_elongation.mean():.3f} | {g.pred_elongation.mean():.3f} | {g.gt_d_eq.mean():.2f} | "
                     f"{g.pred_d_eq.mean():.2f} |")
    L += ["", "## Scoring modes (pooled F1@0.50, TP / FP / FN)", "",
          "| scope | method | " + " | ".join(MODES) + " |", "|---|---|" + "---|" * len(MODES)]
    for scope, d in scopes(modes):
        for method, g in d.groupby("method"):
            cells = []
            for mode in MODES:
                h = g[g["mode"] == mode]
                tp, fp, fn = (int(h[f"{k}_iou050"].sum()) for k in ("tp", "fp", "fn"))
                cells.append(f"{2 * tp / (2 * tp + fp + fn):.3f} ({tp}/{fp}/{fn})")
            L.append(f"| {scope} | {method} | " + " | ".join(cells) + " |")
    L += ["", "## Contacts (pooled over crops): recall / precision, n reference contacts", "",
          "| scope | method | reference | n ref | recall | precision |", "|---|---|---|---|---|---|"]
    for scope, d in scopes(contacts):
        for (method, ref), g in d.groupby(["method", "reference"]):
            tp_p = g["tp_precision"].sum() if "tp_precision" in g else np.nan
            prec = (tp_p / g.n_pred_pairs.sum()) if method != "GT (6-adj)" else g.tp.sum() / g.n_pred_pairs.sum()
            L.append(f"| {scope} | {method} | {ref} | {int(g.n_ref.sum())} | {g.tp.sum() / g.n_ref.sum():.3f} | "
                     f"{prec:.3f} |")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--crops", nargs="*")
    a = ap.parse_args(argv)
    cr = crops()
    if a.crops:
        cr = cr[cr.crop_id.isin(a.crops)]
    res = []
    with ProcessPoolExecutor(a.workers) as ex:
        for c, r in zip(cr.crop_id, ex.map(one_crop, cr.to_dict("records"))):
            print(c, flush=True)
            res.append(r)
    pairs = pd.concat([r["pairs"] for r in res], ignore_index=True)
    modes = pd.concat([r["modes"] for r in res], ignore_index=True)
    contacts = pd.concat([r["contacts"] for r in res], ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(OUT / "pairs.csv.gz", index=False)
    modes.to_csv(OUT / "modes.csv", index=False)
    contacts.to_csv(OUT / "contacts.csv", index=False)
    md = tables(pairs, modes, contacts)
    (OUT / "tables.md").write_text(md, encoding="utf-8")
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
