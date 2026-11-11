#!/usr/bin/env python3
"""Hybrid test (Phase 3 inverted): cpsam_v2 instances as markers for the frozen v3 marker-controlled watershed.

Methods compared on the 12 dev + 52 test crops (size filter off everywhere):
  cellpose:cpsam_v2   existing predictions (C:\\BA_TOS_DATA\\runs\\cellpose\\cpsam_v2)
  v3                  frozen BA-TOS v3, BA-TOS_sel (existing frozen runs)
  hybrid              seeds = cpsam_v2 labels eroded by SEED_EROSION voxels (Euclidean, per label, a priori, not tuned),
                      largest 6-connected piece per label; a label that erosion removes keeps its deepest voxel.
                      Watershed = the frozen v3 script src/batos_v3/patched/e11_run_fixed_mask_watershed_v3.py, run
                      unchanged with the frozen v3 flags: flattened openrec image, Otsu mask (p1-p99, closing x1,
                      fill holes, components >= 500) united with the v3 baseline markers and the seeds (as v3 unites
                      it with its ToS and baseline markers), relief = -EDT of the mask, 6-connectivity. No v3 candidate
                      stage (no ToS markers, no gate) enters the hybrid.
  cpsam_v2_erode1     trivial control: each cpsam_v2 label minus its 6-connected shell (voxels with a 6-neighbour of
                      another label or background; crop faces do not erode).

Steps (batos env):
  regen    The frozen test/scale drivers deleted the v3 intermediates. Re-run frozen v3 (batos_v3.pipeline.run_crop,
           preset v3, freeze hash checked) into C:\\BA_TOS_DATA\\runs\\hybrid\\v3_regen\\<crop>, require the
           BA-TOS_sel segmentation to be array-identical to the frozen one (else fail loudly), then prune to the
           flattened image, baseline markers, ToS markers, segmentation, v3_run.json and commands.log.
  build    seeds + hybrid watershed + erode1 control -> C:\\BA_TOS_DATA\\runs\\hybrid\\{hybrid,erode1}\\<crop>
  evaluate evaluate_v2 no-marker (fraction 0.50) at IoU .25/.50/.75, morphology + boundary (MSD, HD95) on matched
           pairs at IoU 0.50, contacts (GT 6-adjacency and the authors' ContactPoint lists), via the helpers of
           scripts/shape_contacts_cellpose_v3.py. Outputs results/hybrid/{per_crop.csv, pairs.csv.gz, contacts.csv,
           tables.md}; all_runs rows hybrid_<method>/<crop>.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

ROOT = Path(r"C:\BA_TOS_DATA\runs\hybrid")
REGEN = ROOT / "v3_regen"
CP = Path(r"C:\BA_TOS_DATA\runs\cellpose\cpsam_v2")
V3_DEV = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test\runs\v3\x1")
V3_TEST = Path(r"C:\BA_TOS_DATA\runs\test_protocol\v3")
V3_SEG = Path("e23c") / "watershed" / "segmentations" / "batos_selective_distance_c1.tif"
WS_SCRIPT = REPO / "src" / "batos_v3" / "patched" / "e11_run_fixed_mask_watershed_v3.py"
WS_FLAGS = ["--outside-policy", "drop", "--v3-min-component", "500", "--v3-closing-iterations", "1"]
SEED_EROSION = 2  # voxels; fixed a priori (cpsam_v2 boundary sits ~1 voxel outside GT, Phase 2 §5), never tuned
OUT = REPO / "results" / "hybrid"
SIX = ndi.generate_binary_structure(3, 1)
METHODS = ["cellpose:cpsam_v2", "v3", "hybrid", "cpsam_v2_erode1"]


def crops() -> pd.DataFrame:
    from shape_contacts_cellpose_v3 import crops as sc_crops
    c = sc_crops()
    gray = pd.concat([pd.read_csv(REPO / "data" / "metadata" / n)[["crop_id", "gray_path"]]
                      for n in ("crops_manifest.csv", "crops_manifest_test.csv")])
    return c.merge(gray, on="crop_id", how="left", validate="1:1")


def frozen_seg(c) -> Path:
    return (V3_DEV if c["split"] == "dev" else V3_TEST) / c["crop_id"] / V3_SEG


def one(directory: Path, pattern: str) -> Path:
    hits = sorted(directory.glob(pattern))
    if len(hits) != 1:
        raise RuntimeError(f"expected one {pattern} in {directory}, got {len(hits)}")
    return hits[0]


def regen_paths(out: Path) -> dict:
    return {"flat": one(out / "01_openrec", "*_flatten_openrec_r*_8u.tif"),
            "baseline": one(out / "pipe" / "04_markers", "*_batos_marker_components_raw.tif"),
            "tos": out / "e11" / "markers_tos_min" / "markers_uint16.tif",
            "seg": out / V3_SEG}


# ----------------------------------------------------------------------------------------------------------- regen
def _regen_one(c: dict) -> dict:
    from batos_v3.pipeline import StageError, run_crop
    out = REGEN / c["crop_id"]
    rec = {"crop_id": c["crop_id"]}
    if (out / "REGEN_OK.json").exists():
        return json.loads((out / "REGEN_OK.json").read_text(encoding="utf-8"))
    if out.exists():
        shutil.rmtree(out)  # partial run of this driver (our own output)
    t0 = time.time()
    try:
        run_crop(c["crop_id"], Path(c["gray_path"]), out, "v3")
    except StageError as e:
        raise RuntimeError(f"{c['crop_id']}: frozen v3 rerun failed: {e}") from e
    p = regen_paths(out)
    a, b = tifffile.imread(p["seg"]), tifffile.imread(frozen_seg(c))
    identical = a.shape == b.shape and bool(np.array_equal(a, b))
    if not identical:
        raise RuntimeError(f"{c['crop_id']}: regenerated BA-TOS_sel differs from the frozen segmentation")
    keep = {p["flat"], p["baseline"], p["tos"], p["seg"], out / "v3_run.json", out / "commands.log"}
    for f in list(out.rglob("*")):
        if f.is_file() and f not in keep:
            f.unlink()
    rec |= {"identical_to_frozen": identical, "seconds": round(time.time() - t0, 1), "flat": str(p["flat"]),
            "baseline": str(p["baseline"])}
    (out / "REGEN_OK.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
    return rec


def regen(workers: int, only: list[str] | None) -> None:
    from run_test_v3 import check_freeze
    cfg = check_freeze()
    print(f"v3 freeze OK: tag {cfg['git_tag']}", flush=True)
    cr = crops()
    if only:
        cr = cr[cr.crop_id.isin(only)]
    rows = []
    with ProcessPoolExecutor(workers) as ex:
        for r in ex.map(_regen_one, cr.to_dict("records")):
            print(f"{r['crop_id']}: identical={r['identical_to_frozen']} {r.get('seconds')} s", flush=True)
            rows.append(r)
    pd.DataFrame(rows).to_csv(ROOT / "regen.csv", index=False)


# ----------------------------------------------------------------------------------------------------------- build
def shell(lab: np.ndarray) -> np.ndarray:
    """Voxels of a label with a 6-neighbour of a different value (other label or background); faces do not count."""
    s = np.zeros(lab.shape, bool)
    for ax in range(3):
        a = np.moveaxis(lab, ax, 0)
        d = a[1:] != a[:-1]
        sm = np.moveaxis(s, ax, 0)
        sm[1:] |= d
        sm[:-1] |= d
    return s & (lab > 0)


def make_seeds(lab: np.ndarray, k: int) -> tuple[np.ndarray, dict]:
    sh = shell(lab)
    depth = ndi.distance_transform_edt(~sh)  # distance to the nearest shell voxel (0 on the shell)
    depth[lab == 0] = 0
    keep = (lab > 0) & (depth >= k)
    comp, n = ndi.label(keep, SIX)
    seeds = np.zeros(lab.shape, np.int32)
    ids = np.unique(lab[lab > 0])
    if n:
        idx = np.arange(1, n + 1)
        owner = ndi.maximum(lab, comp, idx).astype(np.int64)
        size = ndi.sum(np.ones_like(comp, np.int32), comp, idx)
        d = pd.DataFrame({"comp": idx, "owner": owner, "size": size}).sort_values(["owner", "size"],
                                                                                  ascending=[True, False])
        best = d.drop_duplicates("owner")
        lut = np.zeros(n + 1, np.int32)
        lut[best.comp.to_numpy()] = best.owner.to_numpy()
        seeds = lut[comp]
    present = set(np.unique(seeds[seeds > 0]).tolist())
    lost = [int(i) for i in ids if int(i) not in present]
    for i in lost:  # deepest voxel of a label that erosion removed
        m = lab == i
        z = np.argmax(np.where(m, depth, -1))
        seeds.flat[z] = i
    info = {"n_labels": int(ids.size), "n_lost_to_erosion": len(lost),
            "n_split_labels": int((pd.Series(owner).value_counts() > 1).sum()) if n else 0,
            "seed_voxels": int((seeds > 0).sum())}
    return seeds, info


def _build_one(c: dict) -> dict:
    crop = c["crop_id"]
    lab = tifffile.imread(CP / crop / "labels.tif").astype(np.int32)
    rec = {"crop_id": crop}
    # trivial control
    e_dir = ROOT / "erode1" / crop
    if not (e_dir / "labels.tif").exists():
        e_dir.mkdir(parents=True, exist_ok=True)
        er = np.where(shell(lab), 0, lab).astype(np.uint16 if lab.max() < 65536 else np.uint32)
        tifffile.imwrite(e_dir / "labels.tif", er, compression="zlib")
    # hybrid
    h_dir = ROOT / "hybrid" / crop
    seg = h_dir / "watershed" / "segmentations" / "cpsam_seeds_distance_c1.tif"
    if not seg.exists():
        h_dir.mkdir(parents=True, exist_ok=True)
        seeds, info = make_seeds(lab, SEED_EROSION)
        tifffile.imwrite(h_dir / "seeds.tif", seeds.astype(np.uint16 if seeds.max() < 65536 else np.uint32),
                         compression="zlib")
        p = regen_paths(REGEN / crop)
        cmd = [sys.executable, "-I", str(WS_SCRIPT), "--image", str(p["flat"]), "--mask-marker", str(h_dir / "seeds.tif"),
               "--mask-marker", str(p["baseline"]), *WS_FLAGS, "--marker", f"cpsam_seeds={h_dir / 'seeds.tif'}",
               "--out-dir", str(h_dir / "watershed")]
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        (h_dir / "commands.log").write_text("$ " + " ".join(cmd) + "\n" + r.stdout + r.stderr, encoding="utf-8")
        if r.returncode:
            raise RuntimeError(f"{crop}: watershed failed rc={r.returncode}: {r.stderr[-500:]}")
        info |= {"seconds_watershed": round(time.time() - t0, 1), "seed_erosion": SEED_EROSION, "command": cmd}
        (h_dir / "hybrid_run.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    rec |= json.loads((h_dir / "hybrid_run.json").read_text(encoding="utf-8"))
    rec.pop("command", None)
    return rec


def build(workers: int, only: list[str] | None) -> None:
    cr = crops()
    if only:
        cr = cr[cr.crop_id.isin(only)]
    rows = []
    with ProcessPoolExecutor(workers) as ex:
        for r in ex.map(_build_one, cr.to_dict("records")):
            print(f"{r['crop_id']}: labels {r['n_labels']}, lost {r['n_lost_to_erosion']}, split {r['n_split_labels']}, "
                  f"ws {r['seconds_watershed']} s", flush=True)
            rows.append(r)
    pd.DataFrame(rows).to_csv(ROOT / "build.csv", index=False)


# -------------------------------------------------------------------------------------------------------- evaluate
def pred_path(method: str, c: dict) -> Path:
    return {"cellpose:cpsam_v2": CP / c["crop_id"] / "labels.tif",
            "v3": frozen_seg(c),
            "hybrid": ROOT / "hybrid" / c["crop_id"] / "watershed" / "segmentations" / "cpsam_seeds_distance_c1.tif",
            "cpsam_v2_erode1": ROOT / "erode1" / c["crop_id"] / "labels.tif"}[method]


def _eval_one(c: dict) -> dict:
    from batos_v2.evaluate_v2 import evaluate
    from shape_contacts_cellpose_v3 import adjacency, author_contacts, boundary, valid_labels
    gt = tifffile.imread(c["labels_path"]).astype(np.int64)
    valid = valid_labels(gt)
    ref = {p for p in adjacency(gt) if p[0] in valid and p[1] in valid}
    au = author_contacts(c["scan"])
    au_set = None
    if au is not None:
        sub = au[au.a.isin(valid) & au.b.isin(valid)]
        au_set = {(min(a, b), max(a, b)) for a, b in zip(sub.a, sub.b)}
    base = {k: c[k] for k in ("crop_id", "split", "specimen", "scan", "test_class")}
    rows, pairs, contacts = [], [], []
    for method in METHODS:
        path = pred_path(method, c)
        if not path.exists():
            raise RuntimeError(f"missing prediction {method} {path}")
        pred = tifffile.imread(path).astype(np.int64)
        s, _, morph = evaluate(pred, gt, None, rule="fraction", ignore_overlap_threshold=0.50, morph_iou=0.50)
        rows.append(base | {"method": method, **s})
        pb, gb = ndi.find_objects(pred), ndi.find_objects(gt)
        msd, hd = [], []
        for p, g in zip(morph.prediction_label, morph.gt_label):
            sl = tuple(slice(max(0, min(a.start, b.start) - 2), max(a.stop, b.stop) + 2)
                       for a, b in zip(pb[p - 1], gb[g - 1]))
            m1, m2 = boundary(pred[sl] == p, gt[sl] == g)
            msd.append(m1)
            hd.append(m2)
        pairs.append(morph.assign(msd=msd, hd95=hd, **base, method=method))
        p2g = dict(zip(morph.prediction_label.astype(int), morph.gt_label.astype(int)))
        g2p = {g: p for p, g in p2g.items()}
        padj = adjacency(pred)
        for rname, rset in (("GT 6-adj", ref), ("authors", au_set)):
            if rset is None:
                continue
            rec = sum(1 for a, b in rset if a in g2p and b in g2p and (min(g2p[a], g2p[b]), max(g2p[a], g2p[b])) in padj)
            cand = [(p2g[a], p2g[b]) for a, b in padj if a in p2g and b in p2g]
            prec = sum(1 for a, b in cand if (min(a, b), max(a, b)) in rset)
            contacts.append(base | {"method": method, "reference": rname, "n_ref": len(rset), "tp": rec,
                                    "n_pred_pairs": len(cand), "tp_precision": prec})
    return {"rows": pd.DataFrame(rows), "pairs": pd.concat(pairs, ignore_index=True),
            "contacts": pd.DataFrame(contacts)}


def scopes(df):
    yield "dev (12 EFR)", df[df.split == "dev"]
    t = df[df.split == "test"]
    yield "test: all 52", t
    for s in sorted(t.specimen.unique()):
        yield f"test: {s}", t[t.specimen == s]


def tables(per: pd.DataFrame, pairs: pd.DataFrame, contacts: pd.DataFrame) -> str:
    L = ["# Hybrid test: cpsam_v2 seeds + frozen v3 watershed, vs cpsam_v2, v3 and cpsam_v2 eroded by 1 voxel",
         "", f"No-marker mode (fraction 0.50), size filter off. Seed erosion {SEED_EROSION} voxels (a priori).", "",
         "## Detection (crop-macro mean ± SD (min); pooled TP/FP/FN and F1)", "",
         "| scope | method | macro F1@0.50 | macro F1@0.75 | mIoU@0.50 | pooled TP/FP/FN @0.50 | pooled F1@0.50 | pooled F1@0.75 |",
         "|---|---|---|---|---|---|---|---|"]
    for scope, d in scopes(per):
        for m in METHODS:
            g = d[d.method == m]
            tp, fp, fn = (int(g[f"{k}_iou050"].sum()) for k in ("tp", "fp", "fn"))
            t7, f7, n7 = (int(g[f"{k}_iou075"].sum()) for k in ("tp", "fp", "fn"))
            L.append(f"| {scope} | {m} | {g.f1_iou050.mean():.3f} ± {g.f1_iou050.std(ddof=1):.3f} ({g.f1_iou050.min():.3f}) "
                     f"| {g.f1_iou075.mean():.3f} | {g.mean_matched_iou050.mean():.3f} | {tp}/{fp}/{fn} "
                     f"| {2 * tp / (2 * tp + fp + fn):.3f} | {2 * t7 / (2 * t7 + f7 + n7):.3f} |")
    pairs = pairs.assign(volume_err_pct=100 * pairs.volume_rel_error, surface_err_pct=100 * pairs.surface_area_rel_error)
    key = ["crop_id", "gt_label"]
    nm = pairs.groupby(key).method.nunique()
    common = pairs.set_index(key).loc[nm[nm == len(METHODS)].index].reset_index()
    cols = [("volume_err_pct", "volume err %"), ("d_eq_error", "d_eq err"), ("surface_err_pct", "surface err %"),
            ("sphericity_error", "sphericity err"), ("elongation_error", "elongation err"),
            ("centroid_error", "centroid err"), ("msd", "MSD"), ("hd95", "HD95")]
    L += ["", f"## Shape and boundary on GT grains matched by ALL {len(METHODS)} methods (paired): "
          "mean signed error [median |error|], voxels", "",
          "| scope | method | n | " + " | ".join(n for _, n in cols) + " |", "|---|---|---|" + "---|" * len(cols)]
    for scope, d in scopes(common):
        for m in METHODS:
            g = d[d.method == m]
            cells = []
            for col, _ in cols:
                x = g[col].replace([np.inf, -np.inf], np.nan).dropna()
                cells.append(f"{x.mean():.2f} [{x.median():.2f}]" if col in ("msd", "hd95", "centroid_error")
                             else (f"{x.mean():+.1f} [{x.abs().median():.1f}]" if col.endswith("pct")
                                   else f"{x.mean():+.3f} [{x.abs().median():.3f}]"))
            L.append(f"| {scope} | {m} | {len(g)} | " + " | ".join(cells) + " |")
    L += ["", "## Contacts (pooled): recall / precision", "",
          "| scope | method | vs GT 6-adj (n ref) | vs authors' ContactPoint (n ref) |", "|---|---|---|---|"]
    for scope, d in scopes(contacts):
        for m in METHODS:
            cells = []
            for ref in ("GT 6-adj", "authors"):
                g = d[(d.method == m) & (d.reference == ref)]
                cells.append(f"{g.tp.sum() / g.n_ref.sum():.3f} / {g.tp_precision.sum() / g.n_pred_pairs.sum():.3f} "
                             f"({int(g.n_ref.sum())})" if len(g) else "—")
            L.append(f"| {scope} | {m} | " + " | ".join(cells) + " |")
    return "\n".join(L) + "\n"


def evaluate_all(workers: int) -> None:
    cr = crops()
    res = []
    with ProcessPoolExecutor(workers) as ex:
        for c, r in zip(cr.crop_id, ex.map(_eval_one, cr.to_dict("records"))):
            print(c, flush=True)
            res.append(r)
    per = pd.concat([r["rows"] for r in res], ignore_index=True)
    pairs = pd.concat([r["pairs"] for r in res], ignore_index=True)
    contacts = pd.concat([r["contacts"] for r in res], ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)
    per.to_csv(OUT / "per_crop.csv", index=False)
    pairs.to_csv(OUT / "pairs.csv.gz", index=False)
    contacts.to_csv(OUT / "contacts.csv", index=False)
    md = tables(per, pairs, contacts)
    (OUT / "tables.md").write_text(md, encoding="utf-8")

    from batos_v2.extinction_run import upsert
    from batos_v2.run_crop import git_state
    git = git_state()["commit"]
    rows = []
    for r in per[per.method.isin(["hybrid", "cpsam_v2_erode1"])].itertuples():
        run_json = (ROOT / ("hybrid" if r.method == "hybrid" else "erode1") / r.crop_id /
                    ("hybrid_run.json" if r.method == "hybrid" else "labels.tif"))
        rows.append({"experiment_id": f"hybrid_{r.method}/{r.crop_id}", "method": r.method,
                     "model": "cpsam_v2 + v3 watershed" if r.method == "hybrid" else "cpsam_v2 eroded 1 voxel",
                     "volume": r.scan, "crop": r.crop_id, "device": "cpu (watershed) / cuda (cpsam_v2)",
                     "parameters": f"size_filter=none;eval=evaluate_v2:no_markers_fraction050;split={r.split};"
                                   f"seed_erosion={SEED_EROSION if r.method == 'hybrid' else 'n/a'}",
                     "tp": r.tp_iou050, "fp": r.fp_iou050, "fn": r.fn_iou050, "precision": r.precision_iou050,
                     "recall": r.recall_iou050, "f1_025": r.f1_iou025, "f1_050": r.f1_iou050, "f1_075": r.f1_iou075,
                     "ap50": r.ap_iou050, "ap75": r.ap_iou075, "mean_matched_iou": r.mean_matched_iou050,
                     "n_pred": r.prediction_count_evaluable, "n_gt": r.valid_gt_count, "git_commit": git,
                     "run_json": str(run_json)})
    upsert(REPO / "results" / "all_runs.csv", rows, ["experiment_id", "method", "parameters"])
    print(md)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["regen", "build", "evaluate"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--crops", nargs="*")
    a = ap.parse_args(argv)
    if a.step == "regen":
        regen(a.workers, a.crops)
    elif a.step == "build":
        build(a.workers, a.crops)
    else:
        evaluate_all(a.workers)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
