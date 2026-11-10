#!/usr/bin/env python3
"""Phase 2 robustness: degrade the input and run cpsam_v2 and frozen BA-TOS v3 (no parameter change).

Crops: EFRGP01_00_crop04 (dev) and EORGP01_00_test04 (EO test crop on which v3 works best, F1 0.985; every EO
interior crop is already a v3 capacity failure, 0.40-0.70, so degradation would not be visible there).
Conditions (GT-free, seed 20261109):
  clean
  noise2 / noise5 / noise10 : additive Gaussian noise, sigma = 2 / 5 / 10 % of the grain-air contrast,
                              contrast = difference of the two Otsu class means of the crop (GT-free); clipped to
                              the uint16 range.
  down05                    : 0.5x downsampling by 2x2x2 block mean (100^3). Predictions are upsampled back to
                              200^3 by nearest neighbour (np.repeat) and scored against the original GT.
Steps (batos env):
  prepare  -> C:\\BA_TOS_DATA\\runs\\robustness\\inputs\\<crop>__<cond>.tif + manifest.csv
  v3       -> frozen batos_v3.pipeline.run_crop(preset v3) per input (freeze hash check first), --workers N
  evaluate -> evaluate_v2 no-marker (fraction 0.50) and strict (nothing ignored), filter off
cellpose (cellpose env, between prepare and evaluate):
  python -I scripts\\run_cellpose.py --config configs\\cellpose_cpsam_v2.yaml --tag robustness_cpsam_v2 \\
         --manifest C:\\BA_TOS_DATA\\runs\\robustness\\inputs\\manifest.csv --set extra
Outputs: results/cellpose/robustness/{robustness.csv, tables.md}.
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from skimage.filters import threshold_otsu

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

ROOT = Path(r"C:\BA_TOS_DATA\runs\robustness")
INPUTS = ROOT / "inputs"
CP_ROOT = Path(r"C:\BA_TOS_DATA\runs\cellpose\robustness_cpsam_v2")
OUT = REPO / "results" / "cellpose" / "robustness"
CROPS = ["EFRGP01_00_crop04", "EORGP01_00_test04"]
NOISE = {"noise2": 0.02, "noise5": 0.05, "noise10": 0.10}
CONDS = ["clean", *NOISE, "down05"]
SEED = 20261109
MODES = {"no_markers_fraction050": dict(rule="fraction", ignore_overlap_threshold=0.50),
         "strict_noignore": dict(rule="fraction", ignore_overlap_threshold=1.01)}


def manifest_rows() -> pd.DataFrame:
    rows = []
    for name in ("crops_manifest.csv", "crops_manifest_test.csv"):
        d = pd.read_csv(REPO / "data" / "metadata" / name)
        rows.append(d[d.crop_id.isin(CROPS)][["crop_id", "gray_path", "labels_path"]])
    m = pd.concat(rows).set_index("crop_id")
    if len(m) != len(CROPS):
        raise SystemExit(f"crops missing from manifests: {set(CROPS) - set(m.index)}")
    return m


def prepare() -> None:
    INPUTS.mkdir(parents=True, exist_ok=True)
    m = manifest_rows()
    rows, info = [], []
    for k, crop in enumerate(CROPS):
        g = tifffile.imread(m.loc[crop, "gray_path"])
        x = g.astype(np.float64)
        t = threshold_otsu(x)
        contrast = float(x[x > t].mean() - x[x <= t].mean())
        rng = np.random.default_rng(SEED + k)
        for cond in CONDS:
            if cond == "clean":
                y = g
            elif cond in NOISE:
                sigma = NOISE[cond] * contrast
                y = np.clip(np.rint(x + rng.normal(0.0, sigma, x.shape)), 0, 65535).astype(np.uint16)
            else:
                z, yy, xx = (s // 2 for s in g.shape)
                y = np.rint(x[:2 * z, :2 * yy, :2 * xx].reshape(z, 2, yy, 2, xx, 2).mean(axis=(1, 3, 5))).astype(np.uint16)
            name = f"{crop}__{cond}"
            p = INPUTS / f"{name}.tif"
            tifffile.imwrite(p, y, compression="zlib")
            rows.append({"crop_id": name, "gray_path": str(p)})
            info.append({"crop": crop, "condition": cond, "otsu": float(t), "contrast": contrast,
                         "sigma": NOISE.get(cond, 0.0) * contrast, "shape": "x".join(map(str, y.shape))})
    pd.DataFrame(rows).to_csv(INPUTS / "manifest.csv", index=False)
    pd.DataFrame(info).to_csv(INPUTS / "conditions.csv", index=False)
    print(pd.DataFrame(info).to_string())


def _v3_one(name: str) -> str:
    from batos_v3.pipeline import StageError, run_crop
    out = ROOT / "v3" / name
    if (out / "DONE").exists() or (out / "FAILED.txt").exists():
        return f"{name}: already done"
    try:
        run_crop(name, INPUTS / f"{name}.tif", out, "v3")
        (out / "DONE").write_text("ok", encoding="utf-8")
        return f"{name}: done"
    except StageError as e:
        out.mkdir(parents=True, exist_ok=True)
        (out / "FAILED.txt").write_text(str(e), encoding="utf-8")
        return f"{name}: FAILED {e}"


def run_v3(workers: int) -> None:
    from run_test_v3 import check_freeze
    cfg = check_freeze()
    print(f"v3 freeze OK: tag {cfg['git_tag']}", flush=True)
    names = [f"{c}__{k}" for c in CROPS for k in CONDS]
    with ProcessPoolExecutor(workers) as ex:
        for msg in ex.map(_v3_one, names):
            print(msg, flush=True)


def upsample(pred: np.ndarray, shape) -> np.ndarray:
    if pred.shape == tuple(shape):
        return pred
    up = pred.repeat(2, 0).repeat(2, 1).repeat(2, 2)
    if up.shape != tuple(shape):
        raise SystemExit(f"cannot upsample {pred.shape} to {shape}")
    return up


def evaluate_all() -> None:
    from batos_v2.evaluate_v2 import evaluate
    from batos_v3.pipeline import outputs
    m = manifest_rows()
    rows = []
    for crop in CROPS:
        gt = tifffile.imread(m.loc[crop, "labels_path"])
        for cond in CONDS:
            name = f"{crop}__{cond}"
            preds = {}
            cp = CP_ROOT / name
            preds["cellpose:cpsam_v2"] = (tifffile.imread(cp / "labels.tif") if (cp / "labels.tif").exists()
                                          else None, "FAILED" if (cp / "FAILED.json").exists() else "missing")
            v = ROOT / "v3" / name
            if (v / "DONE").exists():
                shape = tifffile.imread(INPUTS / f"{name}.tif").shape
                preds["v3"] = (outputs(v, shape)["BA-TOS_sel"][0], "ok")
            else:
                preds["v3"] = (None, "FAILED" if (v / "FAILED.txt").exists() else "missing")
            for method, (pred, status) in preds.items():
                if pred is None:
                    if status == "missing":
                        raise SystemExit(f"missing prediction {method} {name}")
                    pred = np.zeros(gt.shape, np.uint16)  # failure scored as an empty segmentation
                else:
                    status = "ok"
                pred = upsample(pred, gt.shape)
                for mode, kw in MODES.items():
                    s, _, _ = evaluate(pred, gt, None, morph_iou=None, **kw)
                    rows.append({"crop": crop, "condition": cond, "method": method, "status": status, "mode": mode,
                                 **{k: s[k] for k in ("f1_iou050", "f1_iou075", "precision_iou050", "recall_iou050",
                                                      "tp_iou050", "fp_iou050", "fn_iou050", "mean_matched_iou050",
                                                      "n_pred_total", "valid_gt_count")}})
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "robustness.csv", index=False)
    L = ["# Robustness: degraded inputs, cpsam_v2 vs BA-TOS v3 (frozen), size filter off", ""]
    for mode in MODES:
        L += [f"## {mode}: F1@0.50 (F1@0.75)", "", "| crop | method | " + " | ".join(CONDS) + " |",
              "|---|---|" + "---|" * len(CONDS)]
        for (crop, method), g in df[df["mode"] == mode].groupby(["crop", "method"]):
            g = g.set_index("condition")
            cells = [f"{g.loc[c, 'f1_iou050']:.3f} ({g.loc[c, 'f1_iou075']:.3f})"
                     + ("" if g.loc[c, "status"] == "ok" else " FAILED") for c in CONDS]
            L.append(f"| {crop} | {method} | " + " | ".join(cells) + " |")
        L.append("")
    (OUT / "tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prepare", "v3", "evaluate"])
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args(argv)
    {"prepare": prepare, "v3": lambda: run_v3(a.workers), "evaluate": evaluate_all}[a.step]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
