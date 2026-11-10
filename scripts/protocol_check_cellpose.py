#!/usr/bin/env python3
"""Phase 2 protocol check: Cellpose vs BA-TOS v3 under the no-marker rule and under a STRICT rule.

no_markers_fraction050 : evaluate_v2 rule=fraction, threshold 0.50 (a prediction whose own voxels lie >= 50 % on
                         excluded GT is ignored) -- the mode used in results/cellpose/tables.md.
strict_noignore        : evaluate_v2 rule=fraction, threshold 1.01 (> 1, so NO prediction is ever ignored: every
                         prediction is matched or counted as FP). The reference protocol is unchanged: valid GT =
                         >= 1000 voxels and not touching a face; excluded GT is not counted as FN; IoU is computed on
                         the domain outside excluded GT (as in all modes).

Predictions: cellpose <runs_root>/<tag>/<crop>/labels.tif; v3 BA-TOS_sel segmentation
(e23c/watershed/segmentations/batos_selective_distance_c1.tif) from the frozen v3 runs (dev: v3_scale_test x1,
test: test_protocol/v3). Size filter off. Writes results/cellpose/protocol_check/{per_crop.csv, tables.md}.
Batos env: C:\\BA_TOS_DATA\\envs\\batos-py312\\Scripts\\python.exe -I scripts\\protocol_check_cellpose.py --tags cpsam_v2
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from batos_v2.evaluate_v2 import evaluate  # noqa: E402
from evaluate_cellpose import RUNS, V3_DEV_RUNS, V3_TEST_RUNS, crops, scopes  # noqa: E402

V3_SEG = Path("e23c") / "watershed" / "segmentations" / "batos_selective_distance_c1.tif"
MODES = {"no_markers_fraction050": dict(rule="fraction", ignore_overlap_threshold=0.50),
         "strict_noignore": dict(rule="fraction", ignore_overlap_threshold=1.01)}
OUT = REPO / "results" / "cellpose" / "protocol_check"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", nargs="+", required=True)
    ap.add_argument("--crops", nargs="*", help="restrict to these crop ids")
    args = ap.parse_args(argv)
    cr = crops()
    if args.crops:
        cr = cr[cr.crop_id.isin(args.crops)]
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "per_crop.csv"
    old = pd.read_csv(cache) if cache.exists() else pd.DataFrame(columns=["crop_id", "method"])
    done = set(zip(old.crop_id, old.method))
    rows = []
    for c in cr.itertuples():
        gt = None
        preds = {f"cellpose:{t}": RUNS / t / c.crop_id / "labels.tif" for t in args.tags}
        preds["v3"] = (V3_DEV_RUNS if c.split == "dev" else V3_TEST_RUNS) / c.crop_id / V3_SEG
        for method, path in preds.items():
            if (c.crop_id, method) in done:
                continue
            if not path.exists():
                raise SystemExit(f"missing prediction {path}")
            gt = tifffile.imread(c.labels_path) if gt is None else gt
            pred = tifffile.imread(path)
            for mode, kw in MODES.items():
                m, _, _ = evaluate(pred, gt, None, morph_iou=None, **kw)
                rows.append({"crop_id": c.crop_id, "split": c.split, "specimen": c.specimen, "scan": c.scan,
                             "test_class": c.test_class, "method": method, "mode": mode, "status": "ok", **m})
        print(c.crop_id, flush=True)
    per = pd.concat([old, pd.DataFrame(rows)], ignore_index=True) if rows else old
    per.to_csv(cache, index=False)

    L = ["# Protocol check: no-marker rule vs STRICT (nothing ignored), size filter off", ""]
    L += ["| scope | method | mode | n | F1@0.50 mean ± SD | min | pooled TP | FP | FN | pooled F1@0.50 | F1@0.75 pooled "
          "| pred total | ignored | pred/GT |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for scope, d in scopes(per):
        for (method, mode), g in d.groupby(["method", "mode"], sort=True):
            tp, fp, fn = (int(g[f"{k}_iou050"].sum()) for k in ("tp", "fp", "fn"))
            tp75, fp75, fn75 = (int(g[f"{k}_iou075"].sum()) for k in ("tp", "fp", "fn"))
            L.append(f"| {scope} | {method} | {mode} | {len(g)} | {g.f1_iou050.mean():.3f} ± "
                     f"{g.f1_iou050.std(ddof=1) if len(g) > 1 else 0:.3f} | {g.f1_iou050.min():.3f} | {tp} | {fp} | {fn} "
                     f"| {2 * tp / (2 * tp + fp + fn):.3f} | {2 * tp75 / (2 * tp75 + fp75 + fn75):.3f} "
                     f"| {int(g.n_pred_total.sum())} | {int(g.n_pred_ignored.sum())} "
                     f"| {g.prediction_count_evaluable.sum() / g.valid_gt_count.sum():.3f} |")
    (OUT / "tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
