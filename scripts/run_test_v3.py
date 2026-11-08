#!/usr/bin/env python3
"""Phase 5e step 4: frozen BA-TOS v3 (configs/batos_paper_v3_scalefree.yaml, tag batos-paper-v3-scalefree) run ONCE
on the 52 held-out test crops (data/metadata/crops_test.json), evaluated exactly like Phase 1f.

frozen and fixflat are NOT re-run: their Phase-1f test outputs (C:\\BA_TOS_DATA\\runs\\test_protocol\\{frozen,fixflat},
eval\\<crop>.{csv,json}; tags batos-paper-v2-frozen / batos-paper-v2-fixflat) are deterministic frozen-code results and
are reused side by side (scripts/aggregate_test_v3.py).

Per crop, sequentially: batos_v3.pipeline preset v3 (GT-free; r_granulo scale; EMPTY = recorded empty result, F1 0;
any other stage failure -> FAILED.json, scored as an empty segmentation, no retry) -> evaluate_v2 on routes ToS / sel x
size filter none / 500 / own (500 s^3) x modes markers_dominant (seed-ignore, primary) / no_markers_fraction050,
IoU .25/.50/.75, valid GT >= 1000 voxels not touching a face -> prune intermediate TIFFs (keep segmentations, markers,
the scale JSON and stage CSV/JSON/logs for forensics).
Outputs: C:\\BA_TOS_DATA\\runs\\test_protocol\\v3\\<crop>, ...\\eval_v3\\<crop>.{csv,json}; log ...\\test_v3.log.
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
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.evaluate_v2 import evaluate  # noqa: E402
from batos_v3.pipeline import StageError, outputs, run_crop  # noqa: E402

META = REPO / "data" / "metadata" / "crops_test.json"
CONFIG = REPO / "configs" / "batos_paper_v3_scalefree.yaml"
ROOT = Path(r"C:\BA_TOS_DATA\runs\test_protocol")
V3_ROOT = ROOT / "v3"
EVAL = ROOT / "eval_v3"
MODES = {"markers_dominant": dict(rule="dominant"),
         "no_markers_fraction050": dict(rule="fraction", ignore_overlap_threshold=0.50)}
KEEP = ("segmentations/tos_e11_min_distance_c1.tif", "segmentations/batos_selective_distance_c1.tif",
        "normalized_markers/tos_e11_min.tif", "normalized_markers/batos_selective.tif")


def log(msg: str) -> None:
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def check_freeze() -> dict:
    """Fail loudly unless every file recorded in the frozen config still has its frozen SHA-256."""
    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    if cfg.get("status") != "frozen_for_test":
        raise RuntimeError(f"{CONFIG} is not frozen_for_test")
    bad = [k for k, v in cfg["code_sha256_at_freeze"].items() if sha256(REPO / k) != v]
    bad += [k for k, v in cfg["recovered_scripts_sha256"].items() if sha256(Path(v["path"])) != v["sha256"]]
    if bad:
        raise RuntimeError(f"code changed since the v3 freeze: {bad}")
    return cfg


def size_filter(pred: np.ndarray, min_size: int) -> np.ndarray:
    ids, counts = np.unique(pred, return_counts=True)
    keep = ids[(ids > 0) & (counts >= min_size)]
    return np.where(np.isin(pred, keep), pred, 0).astype(pred.dtype)


def run_v3(c: dict) -> None:
    out = V3_ROOT / c["crop_id"]
    if any((out / n).exists() for n in ("DONE", "EMPTY", "FAILED.json")):
        return
    if out.exists():
        shutil.rmtree(out)  # partial run of this driver (our own output)
    t0 = time.time()
    try:
        rec = run_crop(c["crop_id"], Path(c["gray_crop_path"]), out, "v3")
        log(f"{c['crop_id']} v3 done in {time.time() - t0:.0f}s | r_granulo {rec.get('r_granulo')} s {rec.get('s')} "
            f"openrec {rec['openrec_radius']} ({rec.get('openrec_implementation')}) thr03 "
            f"{rec.get('stage03_dark_threshold')} gate {rec.get('gate')} markers {rec.get('tos_markers')} "
            f"{'EMPTY' if rec.get('empty_result') else ''}")
    except StageError as e:  # method failure: recorded, evaluated as an empty segmentation (no tuning, no retry)
        out.mkdir(parents=True, exist_ok=True)
        (out / "FAILED.json").write_text(json.dumps({"method": "v3", "error": str(e),
                                                     "seconds": round(time.time() - t0, 1)}, indent=1),
                                         encoding="utf-8")
        log(f"{c['crop_id']} v3 FAILED: {e}")


def evaluate_crop(c: dict) -> None:
    crop = c["crop_id"]
    done = EVAL / f"{crop}.json"
    if done.exists():
        return
    t0 = time.time()
    out = V3_ROOT / crop
    gt = tifffile.imread(c["gt_path"])
    failed = (out / "FAILED.json").exists()
    rec = json.loads((out / ("FAILED.json" if failed else "v3_run.json")).read_text(encoding="utf-8"))
    res = ({r: (np.zeros(gt.shape, np.uint32),) * 2 for r in ("BA-TOS_ToS", "BA-TOS_sel")} if failed
           else outputs(out, gt.shape))
    own = int(rec.get("size_filter") or 500)
    rows = []
    for route, (pred, markers) in res.items():
        for flt, ms in (("none", None), ("500", 500), ("own", own)):
            p = pred if ms is None else size_filter(pred, ms)
            for mode, kw in MODES.items():
                m, _, _ = evaluate(p, gt, markers if mode.startswith("markers") else None, morph_iou=None, **kw)
                rows.append({"crop_id": crop, "specimen": c["specimen"], "scan": c["scan"],
                             "test_class": c["test_class"], "method": "v3", "route": route, "size_filter": flt,
                             "size_filter_voxels": ms, "mode": mode, **m})
    g = rec.get("gate") or {}
    info = {"crop_id": crop, "specimen": c["specimen"], "scan": c["scan"], "test_class": c["test_class"],
            "v3_r_granulo": rec.get("r_granulo"), "v3_r0": rec.get("r0"),
            "v3_specimen_fraction": rec.get("specimen_fraction"), "v3_s": rec.get("s"),
            "v3_openrec_radius": rec.get("openrec_radius"), "v3_openrec_impl": rec.get("openrec_implementation"),
            "v3_stage03_thr": rec.get("stage03_dark_threshold"), "v3_bg_ref": rec.get("bg_ref"),
            "v3_gate_triggered": g.get("gate_triggered"), "v3_mv": g.get("region_voxels_median"),
            "v3_qv": g.get("region_voxels_q75"), "v3_median_limit": g.get("median_limit"),
            "v3_q75_limit": g.get("q75_limit"), "v3_tos_markers": rec.get("tos_markers"),
            "v3_empty": rec.get("empty_result"), "v3_failed": rec.get("error") if failed else None,
            "v3_seconds": rec.get("seconds_total", rec.get("seconds")), "v3_peak_mb": rec.get("peak_memory_mb"),
            "v3_size_filter_own": own, "v3_git": (rec.get("git") or {}).get("commit"),
            "eval_seconds": round(time.time() - t0, 1)}
    EVAL.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(EVAL / f"{crop}.csv", index=False)
    done.write_text(json.dumps(info, indent=1, default=str), encoding="utf-8")
    f = pd.DataFrame(rows)
    sel = f[(f["mode"] == "markers_dominant") & (f.size_filter == "none")]
    log(f"{crop} eval: " + ", ".join(f"v3/{r.route[7:]} {r.f1_iou050:.3f}" for r in sel.itertuples()))


def prune(crop: str) -> None:
    for p in (V3_ROOT / crop).rglob("*.tif"):
        rel = p.relative_to(V3_ROOT / crop).as_posix()
        if not any(rel.endswith(k) for k in KEEP):
            p.unlink()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--crops", nargs="*")
    ap.add_argument("--no-prune", action="store_true")
    a = ap.parse_args(argv)
    cfg = check_freeze()
    log(f"v3 frozen config OK (tag {cfg['git_tag']}, code commit {cfg['git_commit_code']})")
    meta = json.loads(META.read_text(encoding="utf-8"))
    crops = [c for c in meta["crops"] if not a.crops or c["crop_id"] in a.crops]
    V3_ROOT.mkdir(parents=True, exist_ok=True)
    for i, c in enumerate(crops, 1):
        if not Path(c["gt_path"]).exists():
            raise FileNotFoundError(f"test crop missing: {c['gt_path']}")
        log(f"--- {i}/{len(crops)} {c['crop_id']} ({c['test_class']})")
        run_v3(c)
        evaluate_crop(c)
        if not a.no_prune:
            prune(c["crop_id"])
    log("END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
