#!/usr/bin/env python3
"""Phase 5f step 4: frozen / fixflat / v3 / v3.1 on held-out TEST SET 2 (data/metadata/crops_test2.json), run once.

Every method runs through the code that produced test set 1. The test-1 modules are imported unchanged, and only their
path globals are re-pointed to test set 2:
  frozen, fixflat : scripts/run_test_protocol.py (scale_stage e13, run_frozen, run_fixflat, evaluate_crop, prune)
  v3              : scripts/run_test_v3.py (check_freeze = batos-paper-v3-scalefree SHA-256 check, run_v3,
                    evaluate_crop, prune)
  v3.1            : batos_v31.pipeline.run_crop; evaluated by run_test_v3.evaluate_crop (identical protocol), rows
                    relabelled method 'v3.1'. Its freeze check reads configs/batos_paper_v31.yaml.
Per crop, sequentially: e13 scale (per specimen) -> frozen -> fixflat -> eval -> v3 -> eval -> v3.1 -> eval -> prune.
Failures are recorded and scored as empty segmentations (no retry). Completed crops are skipped on restart; a crop's
partial run directory is deleted by the respective run function.
Outputs: C:\\BA_TOS_DATA\\runs\\test2\\{frozen,fixflat,v3,v31,eval,eval_v3,eval_v31,logs}; aggregation:
scripts/aggregate_test2.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
import run_test_protocol as rtp  # noqa: E402
import run_test_v3 as rtv3  # noqa: E402
from batos_v2.run_crop import git_state  # noqa: E402
from batos_v31.pipeline import StageError  # noqa: E402
from batos_v31.pipeline import run_crop as run_v31_crop  # noqa: E402

META = REPO / "data" / "metadata" / "crops_test2.json"
ROOT = Path(r"C:\BA_TOS_DATA\runs\test2")
CONFIG_V31 = REPO / "configs" / "batos_paper_v31.yaml"

rtp.META, rtp.ROOT, rtp.FROZEN_ROOT, rtp.FIX_ROOT = META, ROOT, ROOT / "frozen", ROOT / "fixflat"
rtv3.META = META
log = rtv3.log


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def check_freeze_v31() -> dict:
    cfg = yaml.safe_load(CONFIG_V31.read_text(encoding="utf-8"))
    if cfg.get("status") != "frozen_for_test":
        raise RuntimeError(f"{CONFIG_V31} is not frozen_for_test")
    bad = [k for k, v in cfg["code_sha256_at_freeze"].items() if sha256(REPO / k) != v]
    bad += [k for k, v in cfg["recovered_scripts_sha256"].items() if sha256(Path(v["path"])) != v["sha256"]]
    if bad:
        raise RuntimeError(f"code changed since the v3.1 freeze: {bad}")
    return cfg


def use_v3_paths(method: str) -> None:
    sub = {"v3": ("v3", "eval_v3"), "v3.1": ("v31", "eval_v31")}[method]
    rtv3.V3_ROOT, rtv3.EVAL = ROOT / sub[0], ROOT / sub[1]


def run_v31(c: dict) -> None:
    out = ROOT / "v31" / c["crop_id"]
    if any((out / n).exists() for n in ("DONE", "EMPTY", "FAILED.json")):
        return
    if out.exists():
        shutil.rmtree(out)
    t0 = time.time()
    try:
        rec = run_v31_crop(c["crop_id"], Path(c["gray_crop_path"]), out)
        log(f"{c['crop_id']} v3.1 done in {time.time() - t0:.0f}s | r_granulo {rec.get('r_granulo')} s {rec.get('s')} "
            f"N {rec.get('n_values')} gate {(rec.get('gate') or {}).get('gate_triggered')} markers {rec.get('tos_markers')}"
            f" {'EMPTY' if rec.get('empty_result') else ''}")
    except StageError as e:
        out.mkdir(parents=True, exist_ok=True)
        (out / "FAILED.json").write_text(json.dumps({"method": "v3.1", "error": str(e),
                                                     "seconds": round(time.time() - t0, 1)}, indent=1),
                                         encoding="utf-8")
        log(f"{c['crop_id']} v3.1 FAILED: {e}")


def evaluate_v31(c: dict) -> None:
    use_v3_paths("v3.1")
    csv_path = rtv3.EVAL / f"{c['crop_id']}.csv"
    fresh = not csv_path.exists()
    rtv3.evaluate_crop(c)
    if fresh:
        d = pd.read_csv(csv_path, dtype={"size_filter": str})
        d["method"] = "v3.1"
        d.to_csv(csv_path, index=False)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--crops", nargs="*")
    ap.add_argument("--no-prune", action="store_true")
    a = ap.parse_args(argv)
    cfg3 = rtv3.check_freeze()
    cfg31 = check_freeze_v31()
    log(f"v3 frozen config OK (tag {cfg3['git_tag']}); v3.1 frozen config OK (tag {cfg31['git_tag']}, "
        f"code commit {cfg31['git_commit_code']})")
    meta = json.loads(META.read_text(encoding="utf-8"))
    crops = [c for c in meta["crops"] if not a.crops or c["crop_id"] in a.crops]
    for d in ("frozen", "fixflat", "v3", "v31", "logs", "eval", "eval_v3", "eval_v31"):
        (ROOT / d).mkdir(parents=True, exist_ok=True)
    git = git_state()
    for i, c in enumerate(crops, 1):
        crop = c["crop_id"]
        if not Path(c["gt_path"]).exists():
            raise FileNotFoundError(f"test-2 crop not cut: {c['gt_path']}")
        log(f"--- {i}/{len(crops)} {crop} ({c['test_class']})")
        rtp.scale_stage(meta, c["specimen"])
        rtp.run_frozen(crop)
        rtp.run_fixflat(crop, git)
        rtp.evaluate_crop(crop, c)
        use_v3_paths("v3")
        rtv3.run_v3(c)
        rtv3.evaluate_crop(c)
        run_v31(c)
        evaluate_v31(c)
        if not a.no_prune:
            rtp.prune(crop)
            for m in ("v3", "v3.1"):
                use_v3_paths(m)
                rtv3.prune(crop)
    log("END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
