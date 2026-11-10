#!/usr/bin/env python3
"""Phase 5g check: batos_fullvolume.stream_eval reproduces batos_v2.evaluate_v2.evaluate (markers, rule dominant)
key by key on frozen-v3 dev-crop outputs, with the crop-face reference protocol, small z-slabs (to exercise the
slab accumulation) and size filter none / 500. Fails loudly on any difference."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import tifffile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.evaluate_v2 import evaluate, reference_protocol  # noqa: E402
from batos_fullvolume.stream_eval import accumulate, finalize, seed_values  # noqa: E402
from batos_v3.pipeline import outputs  # noqa: E402

RUNS = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test\runs\v3\x1")
CROPS = Path(r"C:\BA_TOS_DATA\crops")
IDS = ["EFRGP01_00_crop04", "EFRGP01_00_crop05", "EFRGP01_01_crop03"]


def main() -> int:
    n_bad = 0
    for cid in IDS:
        gt = tifffile.imread(CROPS / cid / "labels.tif")
        res = outputs(RUNS / cid, gt.shape)
        valid = reference_protocol(gt)["valid_ids"]
        for route, (pred, mk) in res.items():
            for ms in (None, 500):
                ids, cnt = np.unique(pred, return_counts=True)
                small = ids[(ids > 0) & (cnt < ms)] if ms else np.zeros(0, np.int64)
                p_ref = np.where(np.isin(pred, small), 0, pred) if ms else pred
                ref, _, _ = evaluate(p_ref, gt, mk, rule="dominant", morph_iou=None)
                acc = accumulate(pred, gt, valid, slab=7)
                sp, sg = seed_values(pred, gt, np.argwhere(mk > 0).astype(np.int32))
                got, _, _ = finalize(acc, sp, sg, remove=small)
                diff = {k: (ref[k], got.get(k)) for k in ref
                        if not (ref[k] == got.get(k) or (isinstance(ref[k], float) and np.isnan(ref[k])
                                                          and np.isnan(got.get(k, np.nan)))
                                or (isinstance(ref[k], float) and abs(ref[k] - got.get(k, np.inf)) < 1e-12))}
                tag = f"{cid} {route} filter={ms}"
                print(f"{tag}: F1@.50 ref {ref['f1_iou050']:.4f} stream {got['f1_iou050']:.4f} "
                      f"keys {len(ref)} diff {len(diff)} {diff if diff else ''}", flush=True)
                n_bad += bool(diff)
    print("PASS" if n_bad == 0 else f"FAIL ({n_bad} cases)")
    return int(n_bad > 0)


if __name__ == "__main__":
    raise SystemExit(main())
