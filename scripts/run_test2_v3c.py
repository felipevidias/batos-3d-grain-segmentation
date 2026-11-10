#!/usr/bin/env python3
"""Phase 5i: v3 and v3c on held-out TEST SET 2, run ONCE (only after the v3c freeze, docs/V3C_CANONICAL.md).

Inputs: data/metadata/crops_test2.json (200^3 probes) and data/metadata/crops_test2_v3c.json (v3c windows of native side
round(200 r_probe / 10) at the probe centre), both written by scripts/cut_test2_v3c.py.
Freeze checks before anything runs (fail loudly):
  v3  : configs/batos_paper_v3_scalefree.yaml hashes over LF-normalised bytes, and `git diff` against the tag
        batos-paper-v3-scalefree empty for every recorded file (the raw-byte check of run_test_v3.check_freeze fails on
        a CRLF checkout of identical content);
  v3c : configs/batos_paper_v3c.yaml (LF-normalised hashes, as recorded by scripts/freeze_v3c.py).
Per crop, sequentially:
  v3  on the probe   : run_test_v3.run_v3 (unchanged) -> outputs in ...\\v3\\<crop>
  v3c on the window  : batos_v3c.canonical.run(r_probe from the probe) -> ...\\v3c\\<crop>
  evaluation         : batos_v3c.eval_lean (exactly equal to batos_v2.evaluate_v2), each method against its own
                       window's GT; routes ToS / sel; modes markers_dominant (primary) / no_markers_fraction050;
                       IoU .25/.50/.75; size filters v3 none / 500 / own (500 s^3), v3c none / own (500 canonical);
                       valid-GT rules 'fixed1000' (>= 1000 native voxels; v3 primary) and 'canon1000'
                       (>= round(1000 (r_probe/10)^3) native voxels; v3c primary) both reported for both methods;
  seed statistics    : markers per valid grain, seeded / split (GT = diagnosis only);
  boundary cost      : on GT grains valid ('canon1000') in BOTH windows and matched (IoU > 0.5) by BOTH methods:
                       mean IoU and median signed / absolute relative volume error (global GT label ids are shared);
  prune              : keep segmentations / markers / JSON.
Failures are recorded and scored as empty segmentations (no retry). Completed crops are skipped on restart.
Outputs: C:\\BA_TOS_DATA\\runs\\test2_v3c\\{v3,v3c,eval}; log C:\\BA_TOS_DATA\\runs\\test2_v3c.log.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
import tifffile
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
import run_test_v3 as rtv3  # noqa: E402
from batos_v3.pipeline import outputs as v3_outputs  # noqa: E402
from batos_v3c import canonical, eval_lean  # noqa: E402

META = REPO / "data" / "metadata" / "crops_test2.json"
META_V3C = REPO / "data" / "metadata" / "crops_test2_v3c.json"
CFG_V3 = REPO / "configs" / "batos_paper_v3_scalefree.yaml"
CFG_V3C = REPO / "configs" / "batos_paper_v3c.yaml"
ROOT = Path(r"C:\BA_TOS_DATA\runs\test2_v3c")
MODES = rtv3.MODES
KEEP = rtv3.KEEP
log = rtv3.log


def sha_lf(p: Path) -> str:
    return hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def check_cfg(cfg_path: Path, tag: str | None) -> dict:
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    if cfg.get("status") != "frozen_for_test":
        raise RuntimeError(f"{cfg_path} is not frozen_for_test")
    bad = [k for k, v in cfg["code_sha256_at_freeze"].items() if sha_lf(REPO / k) != v]
    bad += [k for k, v in cfg["recovered_scripts_sha256"].items() if sha_lf(Path(v["path"])) != v["sha256"]]
    if tag:
        d = subprocess.run(["git", "diff", "--name-only", tag, "--", *cfg["code_sha256_at_freeze"]], cwd=REPO,
                           capture_output=True, text=True, check=True).stdout.split()
        bad += [f"git diff vs {tag}: {x}" for x in d]
    if bad:
        raise RuntimeError(f"code changed since the freeze recorded in {cfg_path.name}: {bad}")
    return cfg


def peak_mb() -> float:
    return round(psutil.Process().memory_info().peak_wset / 2 ** 20, 1)


def run_v3c(c: dict, w: dict) -> dict:
    out = ROOT / "v3c" / c["crop_id"]
    for name in ("v3c_run.json", "FAILED.json"):
        if (out / name).exists():
            return json.loads((out / name).read_text(encoding="utf-8"))
    if out.exists():
        shutil.rmtree(out)
    rec = canonical.run(c["crop_id"], Path(w["gray_crop_path"]), out, r_probe=w["r_probe"])
    cr = rec.get("canonical_run") or {}
    log(f"{c['crop_id']} v3c r {rec['r_granulo']:g} native {w['native_side']}^3 canon {rec['canonical_shape']} gate "
        f"{(cr.get('gate') or {}).get('gate_triggered')} markers {cr.get('tos_markers')} {rec.get('seconds_total')}s "
        f"peak {rec.get('peak_memory_mb')} MB {'FAILED ' + rec['error'] if 'error' in rec else ''}")
    return rec


def size_filter(pred: np.ndarray, n: float) -> np.ndarray:
    return canonical.size_filter(pred, int(np.ceil(n)))


def eval_method(method: str, c: dict, gt: np.ndarray, routes: dict, filters: dict, rules: dict) -> tuple[list, dict]:
    rows, protos = [], {k: eval_lean.protocol(gt, v) for k, v in rules.items()}
    for route, (pred, markers) in routes.items():
        for flt, p in filters(pred, route).items():
            for rule, proto in protos.items():
                for mode, kw in MODES.items():
                    m = eval_lean.evaluate(p, gt, markers if mode.startswith("markers") else None,
                                           min_gt_voxels=rules[rule], proto=proto, **kw)
                    rows.append({"crop_id": c["crop_id"], "specimen": c["specimen"], "scan": c["scan"],
                                 "test_class": c["test_class"], "method": method, "route": route, "size_filter": flt,
                                 "gt_rule": rule, "min_gt_voxels": rules[rule], "mode": mode, **m})
    return rows, protos


def evaluate_crop(c: dict, w: dict, v3c_rec: dict) -> None:
    done = ROOT / "eval" / f"{c['crop_id']}.json"
    if done.exists():
        return
    r = float(w["r_probe"])
    rules = {"fixed1000": 1000, "canon1000": int(w["min_gt_voxels"])}
    # v3 on the probe
    out3 = ROOT / "v3" / c["crop_id"]
    gt3 = tifffile.imread(c["gt_path"])
    failed3 = (out3 / "FAILED.json").exists()
    rec3 = json.loads((out3 / ("FAILED.json" if failed3 else "v3_run.json")).read_text(encoding="utf-8"))
    res3 = ({k: (np.zeros(gt3.shape, np.uint32),) * 2 for k in canonical.ROUTES} if failed3
            else v3_outputs(out3, gt3.shape))
    own3 = int(rec3.get("size_filter") or 500)
    rows3, protos3 = eval_method("v3", c, gt3, res3, lambda p, _: {"none": p, "500": size_filter(p, 500),
                                                                     "own": size_filter(p, own3)}, rules)
    # v3c on its window
    out = ROOT / "v3c" / c["crop_id"]
    gt = tifffile.imread(w["gt_path"])
    nat = canonical.outputs(out, gt.shape)  # <= 300^3 on test set 1's r range; all six arrays fit easily
    res = {k: (v[0], v[2]) for k, v in nat.items()}
    rowsc, protosc = eval_method("v3c", c, gt, res, lambda p, route: {"none": p, "own": nat[route][1]}, rules)
    # seeds and boundary cost (sel route)
    seeds = {f"v3_seed_{k}": v for k, v in eval_lean.seed_stats(res3["BA-TOS_sel"][1], gt3, protos3["canon1000"]).items()}
    seeds |= {f"v3c_seed_{k}": v for k, v in eval_lean.seed_stats(res["BA-TOS_sel"][1], gt, protosc["canon1000"]).items()}
    common_valid = np.intersect1d(protos3["canon1000"]["valid_ids"], protosc["canon1000"]["valid_ids"])
    a = eval_lean.matched(res["BA-TOS_sel"][0], gt, common_valid)
    b = eval_lean.matched(res3["BA-TOS_sel"][0], gt3, common_valid)
    common = sorted(set(a) & set(b))
    bnd = {"bnd_n_valid_both_windows": int(len(common_valid)), "bnd_n_matched_v3c": len(a), "bnd_n_matched_v3": len(b),
           "bnd_n_common": len(common)}
    for name, d in (("v3c", a), ("v3", b)):
        iou = np.array([d[k][1] for k in common])
        ve = np.array([d[k][2] for k in common])
        bnd |= {f"bnd_{name}_miou_common": float(iou.mean()) if len(iou) else None,
                f"bnd_{name}_volerr_median": float(np.median(ve)) if len(ve) else None,
                f"bnd_{name}_abs_volerr_median": float(np.median(np.abs(ve))) if len(ve) else None}
    cr = v3c_rec.get("canonical_run") or {}
    g3 = rec3.get("gate") or {}
    info = {"crop_id": c["crop_id"], "specimen": c["specimen"], "scan": c["scan"], "test_class": c["test_class"],
            "r_probe": r, "v3_r_granulo": rec3.get("r_granulo"), "native_side": w["native_side"],
            "shift_zyx": w["shift_zyx"], "min_gt_canon": rules["canon1000"],
            "v3_failed": rec3.get("error") if failed3 else None, "v3_empty": rec3.get("empty_result"),
            "v3_gate": g3.get("gate_triggered"), "v3_tos_markers": rec3.get("tos_markers"), "v3_size_filter_own": own3,
            "v3c_failed": v3c_rec.get("error"), "v3c_gate": (cr.get("gate") or {}).get("gate_triggered"),
            "v3c_tos_markers": cr.get("tos_markers"), "v3c_seconds": v3c_rec.get("seconds_total"),
            "v3c_peak_mb": v3c_rec.get("peak_memory_mb"), **seeds, **bnd, "driver_peak_mb": peak_mb(),
            "evaluator": "batos_v3c.eval_lean (exact)"}
    (ROOT / "eval").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows3 + rowsc).to_csv(ROOT / "eval" / f"{c['crop_id']}.csv", index=False)
    done.write_text(json.dumps(info, indent=1, default=str), encoding="utf-8")
    f = pd.DataFrame(rows3 + rowsc)
    f = f[(f["mode"] == "markers_dominant") & (f.size_filter == "none") & (f.route == "BA-TOS_sel")]
    pr = {(x.method, x.gt_rule): x.f1_iou050 for x in f.itertuples()}
    log(f"{c['crop_id']} eval sel F1@.50 v3 {pr[('v3', 'fixed1000')]:.3f} (canon {pr[('v3', 'canon1000')]:.3f}) | v3c "
        f"{pr[('v3c', 'canon1000')]:.3f} (fixed {pr[('v3c', 'fixed1000')]:.3f}) | r {r:g} S {w['native_side']} | "
        f"mIoU v3c {bnd.get('bnd_v3c_miou_common')} v3 {bnd.get('bnd_v3_miou_common')} | driver {info['driver_peak_mb']} MB")


def prune(crop: str) -> None:
    for base in (ROOT / "v3" / crop, ROOT / "v3c" / crop / "canon"):
        for p in base.rglob("*.tif"):
            if not any(p.as_posix().endswith(k) for k in KEEP):
                p.unlink()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--crops", nargs="*")
    a = ap.parse_args(argv)
    c3 = check_cfg(CFG_V3, "batos-paper-v3-scalefree")
    cc = check_cfg(CFG_V3C, None)
    log(f"v3 frozen config OK (tag {c3['git_tag']}); v3c frozen config OK (tag {cc['git_tag']}, code {cc['git_commit_code']})")
    meta = json.loads(META.read_text(encoding="utf-8"))
    win = {w["crop_id"]: w for w in json.loads(META_V3C.read_text(encoding="utf-8"))["crops"]}
    crops = [c for c in meta["crops"] if not a.crops or c["crop_id"] in a.crops]
    missing = [c["crop_id"] for c in crops if c["crop_id"] not in win or not Path(c["gt_path"]).exists()
               or not Path(win[c["crop_id"]]["gt_path"]).exists()]
    if missing:
        raise FileNotFoundError(f"test-2 probes / v3c windows not cut: {missing}")
    rtv3.V3_ROOT = ROOT / "v3"
    for d in ("v3", "v3c", "eval"):
        (ROOT / d).mkdir(parents=True, exist_ok=True)
    for i, c in enumerate(crops, 1):
        log(f"--- {i}/{len(crops)} {c['crop_id']} ({c['specimen']}, {c['test_class']})")
        rtv3.run_v3(c)
        rec = run_v3c(c, win[c["crop_id"]])
        evaluate_crop(c, win[c["crop_id"]], rec)
        prune(c["crop_id"])
    log("END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
