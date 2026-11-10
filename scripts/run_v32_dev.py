#!/usr/bin/env python3
"""Phase 5h: BA-TOS v3.2 marker-uniqueness hypotheses (H1/H2/H3) on the 12 EFR dev crops and the synthetic scale set
(f = 0.7 / 1.0 / 1.4 / 2.0). Dev + synthetic only; no test crop is read.

Per crop: v3.1 prefix once (batos_v32.pipeline.run_prefix, kept unpruned for later hypotheses), then per hypothesis the
uniqueness rule + v3 tail (run_tail), evaluated exactly as scripts/run_v3_scale_test.py (routes ToS / sel, modes
markers_dominant / no_markers_fraction050, size filter none / own 500 s^3, IoU .25/.50/.75, valid GT >= round(1000 f^3)).
Seed diagnostics per crop (GT used for diagnosis only, after the run): markers per grain, valid grains seeded (>= 1
marker inside), valid grains split (>= 2 markers inside), markers in background.
Identity checks (hypothesis `none`): EFRGP01_00_crop04 (f = 1) vs v3 and EFRGP01_00_crop02_x2 (f = 2) vs v3.1,
SHA-256 of both segmentations.
Outputs: C:\\BA_TOS_DATA\\runs\\v32_dev\\{prefix,runs/<hyp>,eval/<hyp>}; then scripts/aggregate_v32_dev.py.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from batos_v2.evaluate_v2 import evaluate, reference_protocol  # noqa: E402
from batos_v3.pipeline import StageError, outputs  # noqa: E402
from batos_v32.pipeline import PREFIX_DONE, run_prefix, run_tail  # noqa: E402
from run_v3_scale_test import KEEP, MODES, crop_sets, size_filter  # noqa: E402

ROOT = Path(r"C:\BA_TOS_DATA\runs\v32_dev")
V3 = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test\runs\v3")
V31 = Path(r"C:\BA_TOS_DATA\runs\v31_dev\runs")
IDENTITY = [(1.0, "EFRGP01_00_crop04", V3), (2.0, "EFRGP01_00_crop02_x2", V31)]
SEGS = ("e11/watershed/segmentations/tos_e11_min_distance_c1.tif",
        "e23c/watershed/segmentations/batos_selective_distance_c1.tif")


def log(msg: str) -> None:
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def pdir(f: float, crop: str) -> Path:
    return ROOT / "prefix" / f"x{f:g}" / crop


def tdir(hyp: str, f: float, crop: str) -> Path:
    return ROOT / "runs" / hyp / f"x{f:g}" / crop


def ensure_prefix(f: float, c: dict) -> bool:
    out = pdir(f, c["crop_id"])
    if (out / PREFIX_DONE).exists():
        return True
    if (out / "FAILED.json").exists():
        return False
    if out.exists():
        shutil.rmtree(out)
    t0 = time.time()
    try:
        rec = run_prefix(c["crop_id"], Path(c["gray_crop_path"]), out)
        log(f"{c['crop_id']} prefix {time.time() - t0:.0f}s s {rec['s']} peak {rec['peak_memory_mb']:.0f} MB "
            f"{'EMPTY' if rec.get('empty_result') else ''}")
        prune_prefix(out, rec)
        return True
    except StageError as e:
        out.mkdir(parents=True, exist_ok=True)
        (out / "FAILED.json").write_text(json.dumps({"error": str(e)}), encoding="utf-8")
        log(f"{c['crop_id']} prefix FAILED: {e}")
        return False


def prune_prefix(out: Path, rec: dict) -> None:
    keep = {Path(p).resolve() for p in rec["paths"].values()}
    for p in out.rglob("*.tif"):
        if p.resolve() not in keep and p.parent.resolve() != Path(rec["paths"]["tos_dir"]).resolve():
            p.unlink()


def ensure_tail(hyp: str, f: float, c: dict) -> None:
    out = tdir(hyp, f, c["crop_id"])
    if any((out / n).exists() for n in ("DONE", "EMPTY", "FAILED.json")):
        return
    if out.exists():
        shutil.rmtree(out)
    t0 = time.time()
    try:
        rec = run_tail(pdir(f, c["crop_id"]), hyp, out)
        u = rec.get("uniqueness") or {}
        log(f"{c['crop_id']} {hyp} tail {time.time() - t0:.0f}s markers {u.get('markers_in')} -> {u.get('markers_out')} "
            f"gate {rec.get('gate', {}).get('gate_triggered')}")
    except StageError as e:
        out.mkdir(parents=True, exist_ok=True)
        (out / "FAILED.json").write_text(json.dumps({"error": str(e)}), encoding="utf-8")
        log(f"{c['crop_id']} {hyp} tail FAILED: {e}")


def seed_stats(markers: np.ndarray, gt: np.ndarray, min_gt: int) -> dict:
    proto = reference_protocol(gt, min_gt)
    pos = np.flatnonzero(markers.ravel())
    mk, g = markers.ravel()[pos].astype(np.int64), gt.ravel()[pos].astype(np.int64)
    df = pd.DataFrame({"m": mk, "g": g}).value_counts().reset_index(name="n").sort_values(["m", "n"],
                                                                                          ascending=[True, False])
    owner = df.drop_duplicates("m")  # one GT owner per marker label (majority voxel)
    per = owner[owner.g > 0].g.value_counts()
    valid = proto["valid_ids"]
    cnt = per.reindex(valid, fill_value=0).to_numpy()
    n_all = int(len(np.flatnonzero(np.bincount(gt.ravel().astype(np.int64))[1:])))
    return {"n_markers": int(len(owner)), "markers_background": int((owner.g == 0).sum()), "gt_all": n_all,
            "gt_valid": int(len(valid)), "markers_per_grain_all": len(owner) / n_all if n_all else None,
            "markers_per_valid_grain": float(cnt.mean()) if len(cnt) else None,
            "valid_seeded": int((cnt >= 1).sum()), "valid_unseeded": int((cnt == 0).sum()),
            "valid_split": int((cnt >= 2).sum()),
            "frac_seeded": float((cnt >= 1).mean()) if len(cnt) else None,
            "frac_split": float((cnt >= 2).mean()) if len(cnt) else None}


def evaluate_crop(hyp: str, f: float, c: dict) -> None:
    edir = ROOT / "eval" / hyp
    done = edir / f"{c['crop_id']}.json"
    if done.exists():
        return
    out = tdir(hyp, f, c["crop_id"])
    gt = tifffile.imread(c["gt_path"])
    failed = (out / "FAILED.json").exists()
    rec = json.loads((out / ("FAILED.json" if failed else "v3_run.json")).read_text(encoding="utf-8"))
    res = ({r: (np.zeros(gt.shape, np.uint32),) * 2 for r in ("BA-TOS_ToS", "BA-TOS_sel")} if failed
           else outputs(out, gt.shape))
    own = int(rec.get("size_filter") or 500)
    rows, seeds = [], {}
    for route, (pred, markers) in res.items():
        seeds[route] = seed_stats(markers, gt, c["min_gt_voxels"])
        for flt, ms in (("none", None), ("own", own)):
            p = pred if ms is None else size_filter(pred, ms)
            for mode, kw in MODES.items():
                m, _, _ = evaluate(p, gt, markers if mode.startswith("markers") else None, morph_iou=None,
                                   min_gt_voxels=c["min_gt_voxels"], **kw)
                rows.append({"crop_id": c["crop_id"], "dev_crop_id": c["dev_crop_id"], "factor": f,
                             "method": f"v3.2-{hyp}", "route": route, "size_filter": flt, "size_filter_voxels": ms,
                             "mode": mode, **m, **{f"seed_{k}": v for k, v in seeds[route].items()}})
    g = rec.get("gate") or {}
    info = {"crop_id": c["crop_id"], "dev_crop_id": c["dev_crop_id"], "factor": f, "hypothesis": hyp,
            "s": rec.get("s"), "r_granulo": rec.get("r_granulo"), "uniqueness": rec.get("uniqueness"),
            "gate": g.get("gate_triggered"), "mv": g.get("region_voxels_median"), "qv": g.get("region_voxels_q75"),
            "tos_markers_prefix": rec.get("tos_markers_prefix"), "tos_markers": rec.get("tos_markers"),
            "empty": rec.get("empty_result"), "failed": rec.get("error") if failed else None,
            "seconds_total": rec.get("seconds_total"), "peak_mb": rec.get("peak_memory_mb"), "size_filter": own,
            "seed_stats": seeds}
    edir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(edir / f"{c['crop_id']}.csv", index=False)
    done.write_text(json.dumps(info, indent=1, default=str), encoding="utf-8")
    sel = pd.DataFrame(rows)
    sel = sel[(sel["mode"] == "markers_dominant") & (sel.size_filter == "none") & (sel.route == "BA-TOS_sel")].iloc[0]
    st = seeds["BA-TOS_sel"]
    log(f"{c['crop_id']} {hyp} eval sel F1@.50 {sel.f1_iou050:.3f} @.75 {sel.f1_iou075:.3f} | markers/valid grain "
        f"{st['markers_per_valid_grain']:.2f} seeded {st['frac_seeded']:.2f} split {st['frac_split']:.2f}")


def prune_tail(hyp: str, f: float, crop: str) -> None:
    for p in tdir(hyp, f, crop).rglob("*.tif"):
        if not any(p.as_posix().endswith(k) for k in KEEP):
            p.unlink()


def identity_checks(sets: dict) -> None:
    res = {}
    for f, crop, ref in IDENTITY:
        c = next(x for x in sets[f] if x["crop_id"] == crop)
        if not ensure_prefix(f, c):
            raise RuntimeError(f"identity crop {crop}: prefix failed")
        ensure_tail("none", f, c)
        for rel in SEGS:
            a, b = sha(ref / f"x{f:g}" / crop / rel), sha(tdir("none", f, crop) / rel)
            res[f"{crop}/{rel}"] = {"reference": str(ref), "ref_sha256": a, "v3.2-none_sha256": b, "identical": a == b}
        prune_tail("none", f, crop)
    (ROOT / "identity_check.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    ok = all(v["identical"] for v in res.values())
    log(f"identity checks (none vs v3 at f=1, vs v3.1 at f=2): {'BIT-IDENTICAL' if ok else 'DIFFERENT'}")
    if not ok:
        raise RuntimeError("v3.2 prefix/tail split does not reproduce v3 / v3.1")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--hyps", nargs="*", default=["H1"], choices=["none", "H1", "H2", "H3"])
    ap.add_argument("--factors", nargs="*", type=float, default=[1.0, 2.0, 1.4, 0.7])
    ap.add_argument("--crops", nargs="*", help="dev crop ids (default all 12)")
    ap.add_argument("--skip-identity", action="store_true")
    a = ap.parse_args(argv)
    sets = crop_sets()
    ROOT.mkdir(parents=True, exist_ok=True)
    if not a.skip_identity and not (ROOT / "identity_check.json").exists():
        identity_checks(sets)
    for f in a.factors:
        crops = [c for c in sets[f] if not a.crops or c["dev_crop_id"] in a.crops]
        for i, c in enumerate(crops, 1):
            log(f"--- x{f:g} {i}/{len(crops)} {c['crop_id']}")
            ok = ensure_prefix(f, c)
            for h in a.hyps:
                if ok:
                    ensure_tail(h, f, c)
                else:  # prefix crash: scored as an empty segmentation (F1 0), as in the earlier drivers
                    out = tdir(h, f, c["crop_id"])
                    out.mkdir(parents=True, exist_ok=True)
                    shutil.copy(pdir(f, c["crop_id"]) / "FAILED.json", out / "FAILED.json")
                evaluate_crop(h, f, c)
                prune_tail(h, f, c["crop_id"])
    log("END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
