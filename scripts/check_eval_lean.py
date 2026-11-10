#!/usr/bin/env python3
"""Phase 5i: check that batos_v3c.eval_lean reproduces the frozen evaluator exactly before using it.

1. 200^3 crops (the 12 dev crops at f = 1, v3c outputs): eval_lean.evaluate vs batos_v2.evaluate_v2.evaluate(...)[0] for
   both routes, both size filters, both modes; eval_lean.seed_stats vs run_v32_dev.seed_stats; eval_lean.matched vs
   run_v3c_dev.matched (in memory, every key, exact; NaN == NaN).
2. An all-zero prediction (scored-as-empty failure case) on one crop.
3. EFRGP01_00_crop02_x2_canon (400^3): eval_lean rows vs the stored rows written by the frozen evaluator (exact on
   every metric column, from the CSV), and the process peak working set of the lean evaluation.
--x2-only runs part 3 alone in its own process (clean peak). Writes C:\BA_TOS_DATA\runs\v3c_dev\eval_lean_check[_x2].json;
exits 1 on any difference.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
import tifffile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from batos_v2.evaluate_v2 import evaluate as frozen_evaluate, reference_protocol  # noqa: E402
from batos_v3c import canonical, eval_lean  # noqa: E402
from run_v32_dev import seed_stats  # noqa: E402
from run_v3_scale_test import MODES  # noqa: E402
from run_v3c_dev import ROOT, canon_sets, matched, run_dir  # noqa: E402


def same(a, b) -> bool:
    if isinstance(a, float) or isinstance(b, float):
        if a is None or b is None:
            return a is b
        return (math.isnan(a) and math.isnan(b)) or a == b
    return a == b


def diff(a: dict, b: dict) -> list[str]:
    keys = set(a) | set(b)
    return sorted(k for k in keys if k not in a or k not in b or not same(a[k], b[k]))


def main() -> int:
    P = psutil.Process()
    bad, n = [], 0
    sets = canon_sets()
    x2_only = "--x2-only" in sys.argv  # separate process, so the 400^3 peak is not mixed with the frozen 200^3 runs
    for c in ([] if x2_only else sets[1.0]):
        gt = tifffile.imread(c["gt_path"])
        proto = eval_lean.protocol(gt, c["min_gt_voxels"])
        if not np.array_equal(proto["valid_ids"], reference_protocol(gt, c["min_gt_voxels"])["valid_ids"]):
            bad.append(f"{c['crop_id']} valid_ids")
        for route, (pred, own, mk) in canonical.outputs(run_dir(1.0, c["crop_id"]), gt.shape).items():
            for flt, p in (("none", pred), ("own", own)):
                for mode, kw in MODES.items():
                    m = mk if mode.startswith("markers") else None
                    d = diff(frozen_evaluate(p, gt, m, morph_iou=None, min_gt_voxels=c["min_gt_voxels"], **kw)[0],
                             eval_lean.evaluate(p, gt, m, min_gt_voxels=c["min_gt_voxels"], proto=proto, **kw))
                    n += 1
                    if d:
                        bad.append(f"{c['crop_id']} {route} {flt} {mode}: {d}")
            d = diff(seed_stats(mk, gt, c["min_gt_voxels"]), eval_lean.seed_stats(mk, gt, proto))
            if d:
                bad.append(f"{c['crop_id']} {route} seed_stats: {d}")
            if matched(pred, gt, proto["valid_ids"]) != eval_lean.matched(pred, gt, proto["valid_ids"]):
                bad.append(f"{c['crop_id']} {route} matched")
        print(f"{c['crop_id']}: {n} comparisons so far, {len(bad)} differences", flush=True)
    for mode, kw in ([] if x2_only else MODES.items()):
        z = np.zeros(gt.shape, np.uint32)
        n += 1
        d = diff(frozen_evaluate(z, gt, z if mode.startswith("markers") else None, morph_iou=None,
                                 min_gt_voxels=c["min_gt_voxels"], **kw)[0],
                 eval_lean.evaluate(z, gt, z if mode.startswith("markers") else None, min_gt_voxels=c["min_gt_voxels"], **kw))
        if d:
            bad.append(f"empty prediction {mode}: {d}")
    # 400^3: against the rows stored by the frozen evaluator
    c = next(c for c in sets[2.0] if c["crop_id"] == "EFRGP01_00_crop02_x2_canon")
    stored = pd.read_csv(ROOT / "eval" / f"{c['crop_id']}.csv", float_precision="round_trip")
    base = P.memory_info().peak_wset
    gt = tifffile.imread(c["gt_path"])
    proto = eval_lean.protocol(gt, c["min_gt_voxels"])
    for route, (pred, own, mk) in canonical.outputs(run_dir(2.0, c["crop_id"]), gt.shape).items():
        for flt, p in (("none", pred), ("own", own)):
            for mode, kw in MODES.items():
                m = eval_lean.evaluate(p, gt, mk if mode.startswith("markers") else None,
                                       min_gt_voxels=c["min_gt_voxels"], proto=proto, **kw)
                row = stored[(stored.route == route) & (stored.size_filter == flt) & (stored["mode"] == mode)].iloc[0]
                n += 1
                d = [k for k in m if k in row.index and not same(
                    float(m[k]) if isinstance(m[k], (int, float, np.integer, np.floating)) else m[k],
                    float(row[k]) if isinstance(m[k], (int, float, np.integer, np.floating)) else row[k])]
                if d:
                    bad.append(f"{c['crop_id']} {route} {flt} {mode}: {d}")
    peak = P.memory_info().peak_wset / 2 ** 20
    res = {"comparisons": n, "differences": bad, "x2_process_peak_mb": round(peak, 1),
           "process_peak_before_x2_mb": round(base / 2 ** 20, 1)}
    (ROOT / f"eval_lean_check{'_x2' if x2_only else ''}.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
