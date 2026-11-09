#!/usr/bin/env python3
"""Phase 5f step 2: validate BA-TOS v3.1 (count constants floored at v2) on the 12 EFR dev crops and the synthetic
scale set (f = 0.7 / 1.0 / 1.4 / 2.0). Dev + synthetic only; no test crop is read.

* s <= 1 (f = 1.0: r_granulo 10, s = 1; f = 0.7: s = 0.7): v3.1 flags == v3 flags by construction. This is checked per crop
  against the flags recorded in the v3 run (v3_run.json). v3.1 = v3 there, so the v3 evaluations are reused, and one
  dev crop (EFRGP01_00_crop04) is re-run with v3.1 to confirm the segmentations are bit-identical (SHA-256).
* f = 1.4 / 2.0: v3.1 is run (batos_v31.pipeline) and evaluated exactly as in scripts/run_v3_scale_test.py (routes
  ToS / sel, modes markers_dominant / no_markers_fraction050, size filter none / own 500 s^3, IoU .25/.50/.75, valid GT
  >= round(1000 f^3)).
Outputs: C:\\BA_TOS_DATA\\runs\\v31_dev\\{runs,eval}\\...; then scripts/aggregate_v31_dev.py.
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
from batos_v2.evaluate_v2 import evaluate  # noqa: E402
from batos_v31 import constants as C31  # noqa: E402
from batos_v31.pipeline import StageError, outputs, run_crop  # noqa: E402
from run_v3_scale_test import KEEP, MODES, crop_sets, size_filter  # noqa: E402

V3 = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test")
ROOT = Path(r"C:\BA_TOS_DATA\runs\v31_dev")
IDENTITY_CROP = "EFRGP01_00_crop04"


def log(msg: str) -> None:
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def v3_dir(f: float, crop: str) -> Path:
    return V3 / "runs" / "v3" / f"x{f:g}" / crop


def run_dir(f: float, crop: str) -> Path:
    return ROOT / "runs" / f"x{f:g}" / crop


def check_flags(f: float, c: dict) -> dict:
    rec = json.loads((v3_dir(f, c["crop_id"]) / "v3_run.json").read_text(encoding="utf-8"))
    same = C31.stage_flags(rec["s"]) == rec["flags"]
    if f <= 1.0 and not same:
        raise RuntimeError(f"{c['crop_id']}: v3.1 flags differ from v3 at s = {rec['s']}")
    return {"crop_id": c["crop_id"], "factor": f, "s": rec["s"], "flags_identical_to_v3": same}


def ensure_run(f: float, c: dict) -> None:
    out = run_dir(f, c["crop_id"])
    if any((out / n).exists() for n in ("DONE", "EMPTY", "FAILED.json")):
        return
    if out.exists():
        shutil.rmtree(out)
    t0 = time.time()
    try:
        rec = run_crop(c["crop_id"], Path(c["gray_crop_path"]), out)
        log(f"{c['crop_id']} v3.1 done in {time.time() - t0:.0f}s | s {rec['s']} N {rec['n_values']} "
            f"gate {rec.get('gate')} markers {rec.get('tos_markers')} {'EMPTY' if rec.get('empty_result') else ''}")
    except StageError as e:
        out.mkdir(parents=True, exist_ok=True)
        (out / "FAILED.json").write_text(json.dumps({"error": str(e), "seconds": round(time.time() - t0, 1)}),
                                         encoding="utf-8")
        log(f"{c['crop_id']} v3.1 FAILED: {e}")


def evaluate_crop(f: float, c: dict) -> None:
    done = ROOT / "eval" / f"{c['crop_id']}.json"
    if done.exists():
        return
    out = run_dir(f, c["crop_id"])
    gt = tifffile.imread(c["gt_path"])
    failed = (out / "FAILED.json").exists()
    rec = json.loads((out / ("FAILED.json" if failed else "v3_run.json")).read_text(encoding="utf-8"))
    res = ({r: (np.zeros(gt.shape, np.uint32),) * 2 for r in ("BA-TOS_ToS", "BA-TOS_sel")} if failed
           else outputs(out, gt.shape))
    own = int(rec.get("size_filter") or 500)
    rows = []
    for route, (pred, markers) in res.items():
        for flt, ms in (("none", None), ("own", own)):
            p = pred if ms is None else size_filter(pred, ms)
            for mode, kw in MODES.items():
                m, _, _ = evaluate(p, gt, markers if mode.startswith("markers") else None, morph_iou=None,
                                   min_gt_voxels=c["min_gt_voxels"], **kw)
                rows.append({"crop_id": c["crop_id"], "dev_crop_id": c["dev_crop_id"], "factor": f, "method": "v3.1",
                             "route": route, "size_filter": flt, "size_filter_voxels": ms, "mode": mode, **m})
    g = rec.get("gate") or {}
    info = {"crop_id": c["crop_id"], "dev_crop_id": c["dev_crop_id"], "factor": f, "v31_s": rec.get("s"),
            "v31_n_values": rec.get("n_values"), "v31_gate": g.get("gate_triggered"),
            "v31_mv": g.get("region_voxels_median"), "v31_qv": g.get("region_voxels_q75"),
            "v31_tos_markers": rec.get("tos_markers"), "v31_empty": rec.get("empty_result"),
            "v31_failed": rec.get("error") if failed else None, "v31_seconds": rec.get("seconds_total"),
            "v31_peak_mb": rec.get("peak_memory_mb"), "v31_size_filter": own}
    (ROOT / "eval").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(ROOT / "eval" / f"{c['crop_id']}.csv", index=False)
    done.write_text(json.dumps(info, indent=1, default=str), encoding="utf-8")
    sel = pd.DataFrame(rows)
    sel = sel[(sel["mode"] == "markers_dominant") & (sel.size_filter == "none")]
    log(f"{c['crop_id']} eval: " + ", ".join(f"v3.1/{r.route[7:]} {r.f1_iou050:.3f}" for r in sel.itertuples()))


def prune(f: float, crop: str) -> None:
    for p in run_dir(f, crop).rglob("*.tif"):
        if not any(p.as_posix().endswith(k) for k in KEEP):
            p.unlink()


def identity_check(sets: dict) -> None:
    c = next(x for x in sets[1.0] if x["crop_id"] == IDENTITY_CROP)
    ensure_run(1.0, c)
    res = {}
    for rel in ("e11/watershed/segmentations/tos_e11_min_distance_c1.tif",
                "e23c/watershed/segmentations/batos_selective_distance_c1.tif"):
        a, b = sha(v3_dir(1.0, IDENTITY_CROP) / rel), sha(run_dir(1.0, IDENTITY_CROP) / rel)
        res[rel] = {"v3": a, "v3.1": b, "identical": a == b}
    (ROOT / "identity_check.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    ok = all(v["identical"] for v in res.values())
    log(f"identity check {IDENTITY_CROP} (v3 vs v3.1 at s = 1): {'BIT-IDENTICAL' if ok else 'DIFFERENT'}")
    if not ok:
        raise RuntimeError("v3.1 differs from v3 at s = 1")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--factors", nargs="*", type=float, default=[1.4, 2.0])
    ap.add_argument("--crops", nargs="*", help="dev crop ids (default all 12)")
    ap.add_argument("--no-prune", action="store_true")
    a = ap.parse_args(argv)
    sets = crop_sets()
    ROOT.mkdir(parents=True, exist_ok=True)
    flags = [check_flags(f, c) for f in sorted(sets) for c in sets[f]]
    pd.DataFrame(flags).to_csv(ROOT / "flag_identity.csv", index=False)
    log(f"flags identical to v3 on all {sum(1 for x in flags if x['factor'] <= 1.0)} crops with s <= 1; "
        f"differ on {sum(1 for x in flags if not x['flags_identical_to_v3'])} crops with s > 1")
    identity_check(sets)
    for f in a.factors:
        crops = [c for c in sets[f] if not a.crops or c["dev_crop_id"] in a.crops]
        for i, c in enumerate(crops, 1):
            log(f"--- x{f:g} {i}/{len(crops)} {c['crop_id']}")
            ensure_run(f, c)
            evaluate_crop(f, c)
            if not a.no_prune:
                prune(f, c["crop_id"])
    log("END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
