#!/usr/bin/env python3
"""Phase 5i: BA-TOS v3c (canonical-scale normalization) on the 12 EFR dev crops and the synthetic scale set.

Dev + synthetic only; no test crop is read. One crop at a time.
- Identity first: v3c on EFRGP01_00_crop04 (r_granulo 10, f = 1, frozen e13 openrec radius 15) must reproduce the
  existing batos-paper-v2-fixflat segmentations bit-for-bit.
- Then per factor (default 1.0, 2.0, 1.4, 0.7) and crop: run, evaluate (protocol of run_v3_scale_test.evaluate_crop:
  native GT, valid GT >= round(1000 f^3), routes ToS/sel, filter none/own, modes markers_dominant / no_markers_fraction050),
  record seed statistics (GT = diagnosis only) and the boundary cost against v3 on grains matched by both.
- Prune the canonical intermediates; keep the native outputs and the canonical segmentations/markers.
Outputs: C:\\BA_TOS_DATA\\runs\\v3c_dev\\{runs\\x<f>\\<crop>, eval}; log C:\\BA_TOS_DATA\\runs\\v3c_dev.log.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from batos_v2.evaluate_v2 import evaluate  # noqa: E402
from batos_v3c import canonical  # noqa: E402
from run_v32_dev import seed_stats  # noqa: E402
from run_v3_scale_test import MODES, crop_sets  # noqa: E402

ROOT = Path(r"C:\BA_TOS_DATA\runs\v3c_dev")
V3RUNS = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test\runs\v3")
DEV_FIXFLAT = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_fixflat\A")
KEEP = ("segmentations/tos_e11_min_distance_c1.tif", "segmentations/batos_selective_distance_c1.tif",
        "normalized_markers/tos_e11_min.tif", "normalized_markers/batos_selective.tif")


def log(msg: str) -> None:
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def run_dir(f: float, crop: str) -> Path:
    return ROOT / "runs" / f"x{f:g}" / crop


def ensure(f: float, c: dict) -> dict:
    out = run_dir(f, c["crop_id"])
    for name in ("v3c_run.json", "FAILED.json"):
        if (out / name).exists():
            return json.loads((out / name).read_text(encoding="utf-8"))
    if out.exists():
        shutil.rmtree(out)
    rec = canonical.run(c["crop_id"], Path(c["gray_crop_path"]), out)
    cr = rec.get("canonical_run") or {}
    log(f"{c['crop_id']} v3c r {rec['r_granulo']:g} f {rec['f']:.3f} canon {rec['canonical_shape']} "
        f"gate {(cr.get('gate') or {}).get('gate_triggered')} markers {cr.get('tos_markers')} "
        f"{rec.get('seconds_total')}s peak {rec.get('peak_memory_mb')} MB {'FAILED ' + rec['error'] if 'error' in rec else ''}")
    return rec


def identity_check() -> None:
    p = ROOT / "identity_check.json"
    if p.exists():
        return
    c = next(c for c in crop_sets()[1.0] if c["crop_id"] == "EFRGP01_00_crop04")
    rec = ensure(1.0, c)
    out = run_dir(1.0, c["crop_id"])
    res = {"crop_id": c["crop_id"], "r_granulo": rec["r_granulo"], "resampled": rec.get("resampled")}
    for route, (ws, label) in canonical.SEGMENTATIONS.items():
        ref = tifffile.imread(DEV_FIXFLAT / c["crop_id"] / ws / "segmentations" / f"{label}_distance_c1.tif")
        new = tifffile.imread(out / "canon" / "run" / ws / "segmentations" / f"{label}_distance_c1.tif")
        nat = tifffile.imread(out / "native" / f"{route}_labels.tif")
        res[route] = {"canonical_identical": bool(np.array_equal(ref, new)), "native_identical": bool(np.array_equal(ref, nat))}
    p.write_text(json.dumps(res, indent=1), encoding="utf-8")
    log(f"identity: {res}")
    if not all(v["canonical_identical"] and v["native_identical"] for k, v in res.items() if k in canonical.ROUTES):
        raise SystemExit("identity check FAILED: v3c at f = 1 does not reproduce batos-paper-v2-fixflat")


def matched(pred: np.ndarray, gt: np.ndarray, valid: np.ndarray) -> dict[int, tuple[int, float, float]]:
    """valid GT id -> (pred id, IoU, relative volume error) for IoU > 0.5."""
    p, g = pred.ravel().astype(np.int64), gt.ravel().astype(np.int64)
    sp, sg = np.bincount(p), np.bincount(g)
    both = (p > 0) & (g > 0)
    m = int(g.max()) + 1
    u, inter = np.unique(p[both] * m + g[both], return_counts=True)
    pu, gu = u // m, u % m
    iou = inter / (sp[pu] + sg[gu] - inter)
    vs = set(valid.tolist())
    return {int(b): (int(a), float(i), float(sp[a] / sg[b] - 1)) for a, b, i in zip(pu, gu, iou)
            if i > 0.5 and int(b) in vs}


def boundary_cost(v3c: np.ndarray, v3: np.ndarray | None, gt: np.ndarray, valid: np.ndarray) -> dict:
    a = matched(v3c, gt, valid)
    if v3 is None:
        return {"n_matched_v3c": len(a)}
    b = matched(v3, gt, valid)
    common = sorted(set(a) & set(b))
    out = {"n_matched_v3c": len(a), "n_matched_v3": len(b), "n_common": len(common)}
    for name, d in (("v3c", a), ("v3", b)):
        iou = np.array([d[k][1] for k in common])
        ve = np.array([d[k][2] for k in common])
        out |= {f"{name}_miou_common": float(iou.mean()) if len(iou) else None,
                f"{name}_volerr_median": float(np.median(ve)) if len(ve) else None,
                f"{name}_abs_volerr_median": float(np.median(np.abs(ve))) if len(ve) else None}
    return out


def evaluate_crop(f: float, c: dict, rec: dict) -> None:
    done = ROOT / "eval" / f"{c['crop_id']}.json"
    if done.exists():
        return
    gt = tifffile.imread(c["gt_path"])
    out = run_dir(f, c["crop_id"])
    res = canonical.outputs(out, gt.shape)
    rows = []
    for route, (pred, pred_own, markers) in res.items():
        for flt, p in (("none", pred), ("own", pred_own)):
            for mode, kw in MODES.items():
                m, _, _ = evaluate(p, gt, markers if mode.startswith("markers") else None, morph_iou=None,
                                   min_gt_voxels=c["min_gt_voxels"], **kw)
                rows.append({"crop_id": c["crop_id"], "dev_crop_id": c["dev_crop_id"], "factor": f, "method": "v3c",
                             "route": route, "size_filter": flt,
                             "size_filter_voxels": None if flt == "none" else canonical.SIZE_FILTER / rec["f"] ** 3,
                             "mode": mode, **m})
    from batos_v2.evaluate_v2 import reference_protocol
    valid = np.asarray(reference_protocol(gt, c["min_gt_voxels"])["valid_ids"])
    v3p = V3RUNS / f"x{f:g}" / c["crop_id"] / "e23c" / "watershed" / "segmentations" / "batos_selective_distance_c1.tif"
    sel, _, sel_mk = res["BA-TOS_sel"]
    cr = rec.get("canonical_run") or {}
    info = {"crop_id": c["crop_id"], "dev_crop_id": c["dev_crop_id"], "factor": f, "min_gt_voxels": c["min_gt_voxels"],
            "r_granulo": rec.get("r_granulo"), "f_canon": rec.get("f"), "canonical_shape": rec.get("canonical_shape"),
            "resampled": rec.get("resampled"), "failed": rec.get("error"), "gate": (cr.get("gate") or {}).get("gate_triggered"),
            "mv": (cr.get("gate") or {}).get("region_voxels_median"), "qv": (cr.get("gate") or {}).get("region_voxels_q75"),
            "tos_markers": cr.get("tos_markers"), "seconds": rec.get("seconds_total"), "peak_mb": rec.get("peak_memory_mb"),
            **{f"seed_{k}": v for k, v in seed_stats(sel_mk, gt, c["min_gt_voxels"]).items()},
            **{f"bnd_{k}": v for k, v in boundary_cost(sel, tifffile.imread(v3p) if v3p.exists() else None, gt,
                                                        valid).items()}}
    (ROOT / "eval").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(ROOT / "eval" / f"{c['crop_id']}.csv", index=False)
    done.write_text(json.dumps(info, indent=1, default=str), encoding="utf-8")
    r = pd.DataFrame(rows)
    r = r[(r["mode"] == "markers_dominant") & (r.size_filter == "none")].set_index("route").f1_iou050
    log(f"{c['crop_id']} eval sel {r['BA-TOS_sel']:.3f} ToS {r['BA-TOS_ToS']:.3f} | markers/valid grain "
        f"{info['seed_markers_per_valid_grain']} split {info['seed_valid_split']} | mIoU v3c "
        f"{info.get('bnd_v3c_miou_common')} v3 {info.get('bnd_v3_miou_common')}")


def prune(f: float, crop: str) -> None:
    for p in (run_dir(f, crop) / "canon").rglob("*.tif"):
        if not any(p.as_posix().endswith(k) for k in KEEP):
            p.unlink()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--factors", nargs="*", type=float, default=[1.0, 2.0, 1.4, 0.7])
    ap.add_argument("--crops", nargs="*", help="dev crop ids (default all 12)")
    a = ap.parse_args(argv)
    ROOT.mkdir(parents=True, exist_ok=True)
    identity_check()
    sets = crop_sets()
    for f in a.factors:
        crops = [c for c in sets[f] if not a.crops or c["dev_crop_id"] in a.crops]
        for i, c in enumerate(crops, 1):
            log(f"--- x{f:g} {i}/{len(crops)} {c['crop_id']}")
            rec = ensure(f, c)
            evaluate_crop(f, c, rec)
            prune(f, c["crop_id"])
    log("END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
