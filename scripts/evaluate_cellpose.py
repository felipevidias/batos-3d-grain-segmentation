#!/usr/bin/env python3
"""Phase 2: evaluate Cellpose predictions and compare with BA-TOS v3, both in the NO-MARKER mode.

Runs in the batos env:

    C:\\BA_TOS_DATA\\envs\\batos-py312\\Scripts\\python.exe -I scripts\\evaluate_cellpose.py --tags cpsam_v2 [...]

1. Per crop and tag: <runs_root>/<tag>/<crop>/labels.tif (FAILED.json -> empty segmentation, F1 0; a missing
   prediction is an error) -> batos_v2.evaluate_v2 with rule=fraction, ignore_overlap_threshold=0.50 (= the
   no-marker mode used for v3, `no_markers_fraction050`), IoU 0.25/0.50/0.75, size filter none and 500 (labels
   < 500 voxels set to 0, no relabelling, same function as run_test_v3.py). Cached in
   <runs_root>/<tag>/eval/<crop>.csv (delete to recompute).
2. v3 reference rows are READ from the frozen v3 evaluations (no re-evaluation): dev crops from
   C:\\BA_TOS_DATA\\runs\\v3_scale_test\\eval\\<crop>.csv (size filter `own` = 500 at s = 1), test crops from
   C:\\BA_TOS_DATA\\runs\\test_protocol\\eval_v3\\<crop>.csv; method v3, route BA-TOS_sel, mode no_markers_fraction050.
3. Writes results/cellpose/{per_crop.csv, crop_macro.csv, pooled.csv, runtime.csv, tables.md} and upserts rows
   `cellpose_<tag>/<crop>` into results/all_runs.csv. Crop-macro and pooled tables are separate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.evaluate_v2 import evaluate  # noqa: E402
from batos_v2.extinction_run import upsert  # noqa: E402
from batos_v2.run_crop import git_state  # noqa: E402

RUNS = Path(r"C:\BA_TOS_DATA\runs\cellpose")
V3_DEV_EVAL = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test\eval")
V3_DEV_RUNS = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test\runs\v3\x1")
V3_TEST_EVAL = Path(r"C:\BA_TOS_DATA\runs\test_protocol\eval_v3")
V3_TEST_RUNS = Path(r"C:\BA_TOS_DATA\runs\test_protocol\v3")
OUT = REPO / "results" / "cellpose"
MODE = "no_markers_fraction050"
EVAL_KW = dict(rule="fraction", ignore_overlap_threshold=0.50, morph_iou=None)
TS = ("025", "050", "075")


def crops() -> pd.DataFrame:
    dev = pd.read_csv(REPO / "data" / "metadata" / "crops_manifest.csv")
    dev = dev.assign(split="dev", specimen="EFR", test_class="dev")
    test = pd.read_csv(REPO / "data" / "metadata" / "crops_manifest_test.csv").assign(split="test")
    cols = ["crop_id", "split", "specimen", "scan", "test_class", "labels_path"]
    return pd.concat([dev[cols], test[cols]], ignore_index=True)


def size_filter(pred: np.ndarray, min_size: int) -> np.ndarray:
    ids, counts = np.unique(pred, return_counts=True)
    keep = ids[(ids > 0) & (counts >= min_size)]
    return np.where(np.isin(pred, keep), pred, 0).astype(pred.dtype)


def eval_tag(tag: str, c: pd.Series) -> pd.DataFrame:
    cache = RUNS / tag / "eval" / f"{c.crop_id}.csv"
    if cache.exists():
        return pd.read_csv(cache, dtype={"size_filter": str})
    d = RUNS / tag / c.crop_id
    gt = tifffile.imread(c.labels_path)
    if (d / "FAILED.json").exists():
        pred, status = np.zeros(gt.shape, np.uint16), "failed"
    elif (d / "labels.tif").exists() and (d / "run.json").exists():
        pred, status = tifffile.imread(d / "labels.tif"), "ok"
    else:
        raise SystemExit(f"missing prediction for {tag}/{c.crop_id} (no labels.tif+run.json, no FAILED.json)")
    if pred.shape != gt.shape:
        raise SystemExit(f"shape mismatch {tag}/{c.crop_id}: {pred.shape} vs {gt.shape}")
    rows = []
    for flt, ms in (("none", None), ("500", 500)):
        m, _, _ = evaluate(pred if ms is None else size_filter(pred, ms), gt, None, **EVAL_KW)
        rows.append({"crop_id": c.crop_id, "split": c.split, "specimen": c.specimen, "scan": c.scan,
                     "test_class": c.test_class, "method": f"cellpose:{tag}", "status": status, "size_filter": flt,
                     "mode": MODE, **m})
    df = pd.DataFrame(rows)
    cache.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False)
    return df


def v3_rows(c: pd.Series) -> pd.DataFrame:
    path = (V3_DEV_EVAL if c.split == "dev" else V3_TEST_EVAL) / f"{c.crop_id}.csv"
    if not path.exists():
        raise SystemExit(f"missing v3 evaluation {path}")
    d = pd.read_csv(path, dtype={"size_filter": str})
    d = d[(d.method == "v3") & (d.route == "BA-TOS_sel") & (d["mode"] == MODE)]
    if c.split == "dev":
        d = d[d.size_filter.isin(["none", "own"])].replace({"size_filter": {"own": "500"}})
    else:
        d = d[d.size_filter.isin(["none", "500"])]
    if len(d) != 2:
        raise SystemExit(f"expected 2 v3 rows for {c.crop_id}, got {len(d)}")
    d = d.drop(columns=[k for k in ("specimen", "scan", "test_class", "route") if k in d.columns])
    return d.assign(split=c.split, specimen=c.specimen, scan=c.scan, test_class=c.test_class, method="v3",
                    status="ok")


def scopes(df: pd.DataFrame):
    dev, test = df[df.split == "dev"], df[df.split == "test"]
    yield "dev (12 EFR)", dev
    yield "test: all 52", test
    yield "test: EFR_internal", test[test.specimen == "EFR"]
    yield "test: external", test[test.specimen != "EFR"]
    for s in sorted(test.specimen.unique()):
        yield f"test: {s}", test[test.specimen == s]


def aggregate(per: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    macro, pooled = [], []
    for scope, d in scopes(per):
        for (method, flt), g in d.groupby(["method", "size_filter"], sort=False):
            base = {"scope": scope, "method": method, "size_filter": flt, "n_crops": len(g)}
            r = dict(base)
            for t in TS:
                f = g[f"f1_iou{t}"]
                r |= {f"f1_{t}_mean": f.mean(), f"f1_{t}_sd": f.std(ddof=1) if len(f) > 1 else np.nan,
                      f"f1_{t}_median": f.median(), f"f1_{t}_min": f.min()}
            r["min_crop"] = g.loc[g.f1_iou050.idxmin(), "crop_id"]
            r["n_below_090"] = int((g.f1_iou050 < 0.90).sum())
            r["mean_matched_iou_050"] = g.mean_matched_iou050.mean()
            r["n_failed"] = int((g.status == "failed").sum())
            macro.append(r)
            p = dict(base) | {"valid_gt": int(g.valid_gt_count.sum()), "pred_eval": int(g.prediction_count_evaluable.sum())}
            for t in TS:
                tp, fp, fn = (int(g[f"{k}_iou{t}"].sum()) for k in ("tp", "fp", "fn"))
                pr = tp / (tp + fp) if tp + fp else 0.0
                rc = tp / (tp + fn) if tp + fn else 0.0
                p |= {f"tp_{t}": tp, f"fp_{t}": fp, f"fn_{t}": fn, f"precision_{t}": pr, f"recall_{t}": rc,
                      f"f1_{t}": 2 * tp / (2 * tp + fp + fn) if tp else 0.0}
            pooled.append(p)
    return pd.DataFrame(macro), pd.DataFrame(pooled)


def runtime_table(tags: list[str], cr: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in cr.itertuples():
        for tag in tags:
            j = RUNS / tag / c.crop_id / "run.json"
            if j.exists():
                r = json.loads(j.read_text(encoding="utf-8"))
                rows.append({"crop_id": c.crop_id, "split": c.split, "method": f"cellpose:{tag}",
                             "device": f"cuda ({r['environment']['device_name']})", "seconds": r["runtime_seconds"],
                             "peak_vram_reserved_mb": r["peak_vram_reserved_mb"],
                             "peak_vram_allocated_mb": r["peak_vram_allocated_mb"],
                             "peak_host_mb": r["peak_host_rss_mb"], "n_labels": r["n_labels"], "run_json": str(j)})
        j = (V3_DEV_RUNS if c.split == "dev" else V3_TEST_RUNS) / c.crop_id / "v3_run.json"
        if j.exists():
            r = json.loads(j.read_text(encoding="utf-8"))
            rows.append({"crop_id": c.crop_id, "split": c.split, "method": "v3", "device": "cpu (AMD Ryzen 7 5700X)",
                         "seconds": r.get("seconds_total"), "peak_vram_reserved_mb": np.nan,
                         "peak_vram_allocated_mb": np.nan, "peak_host_mb": r.get("peak_memory_mb"),
                         "n_labels": np.nan, "run_json": str(j)})
    return pd.DataFrame(rows)


def tables_md(macro: pd.DataFrame, pooled: pd.DataFrame, rt: pd.DataFrame) -> str:
    L = ["# Phase 2 — Cellpose vs BA-TOS v3, no-marker evaluation (evaluate_v2, fraction 0.50)", "",
         "v3 = BA-TOS_sel route of `batos-paper-v3-scalefree` (gate never fires on these crops, so sel = ToS).",
         "Crop-macro and pooled tables are separate. Size filter 500 = labels < 500 voxels removed, no relabelling.", ""]
    for flt in ("none", "500"):
        L += [f"## Crop-macro F1, size filter {flt}", "",
              "| scope | method | n | F1@0.25 | F1@0.50 mean ± SD | median | min (crop) | # < 0.90 | F1@0.75 | mIoU@0.50 | failed |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in macro[macro.size_filter == flt].itertuples():
            L.append(f"| {r.scope} | {r.method} | {r.n_crops} | {r.f1_025_mean:.3f} | {r.f1_050_mean:.3f} ± "
                     f"{r.f1_050_sd:.3f} | {r.f1_050_median:.3f} | {r.f1_050_min:.3f} ({r.min_crop}) | {r.n_below_090} "
                     f"| {r.f1_075_mean:.3f} | {r.mean_matched_iou_050:.3f} | {r.n_failed} |")
        L += ["", f"## Pooled, size filter {flt}", "",
              "| scope | method | GT | pred | TP@0.50 | FP@0.50 | FN@0.50 | P@0.50 | R@0.50 | F1@0.25 | F1@0.50 | F1@0.75 |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in pooled[pooled.size_filter == flt].itertuples():
            L.append(f"| {r.scope} | {r.method} | {r.valid_gt} | {r.pred_eval} | {r.tp_050} | {r.fp_050} | {r.fn_050} "
                     f"| {r.precision_050:.3f} | {r.recall_050:.3f} | {r.f1_025:.3f} | {r.f1_050:.3f} | {r.f1_075:.3f} |")
        L.append("")
    L += ["## Runtime and memory per crop (200³)", "",
          "| split | method | device | n | seconds median [min–max] | peak VRAM reserved MB (max) | peak host MB (max) |",
          "|---|---|---|---|---|---|---|"]
    for (split, method, device), g in rt.groupby(["split", "method", "device"], sort=False):
        v = g.peak_vram_reserved_mb.max()
        L.append(f"| {split} | {method} | {device} | {len(g)} | {g.seconds.median():.0f} [{g.seconds.min():.0f}–"
                 f"{g.seconds.max():.0f}] | {'—' if pd.isna(v) else f'{v:.0f}'} | {g.peak_host_mb.max():.0f} |")
    L += ["", "Cellpose seconds = model.eval only (model load excluded), GPU. v3 seconds = whole pipeline, CPU, "
          "not contention-controlled (from v3_run.json)."]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tags", nargs="+", required=True, help="cellpose run tags under C:\\BA_TOS_DATA\\runs\\cellpose")
    ap.add_argument("--split", choices=["dev", "test", "both"], default="both")
    ap.add_argument("--out", type=Path, default=OUT, help="output folder for the tables (default results/cellpose)")
    args = ap.parse_args(argv)
    out = args.out
    cr = crops()
    if args.split != "both":
        cr = cr[cr.split == args.split]
    parts = []
    for c in cr.itertuples(index=False):
        c = pd.Series(c._asdict())
        for tag in args.tags:
            parts.append(eval_tag(tag, c))
        parts.append(v3_rows(c))
    per = pd.concat(parts, ignore_index=True)
    macro, pooled = aggregate(per)
    rt = runtime_table(args.tags, cr)
    out.mkdir(parents=True, exist_ok=True)
    per.to_csv(out / "per_crop.csv", index=False)
    macro.to_csv(out / "crop_macro.csv", index=False)
    pooled.to_csv(out / "pooled.csv", index=False)
    rt.to_csv(out / "runtime.csv", index=False)
    (out / "tables.md").write_text(tables_md(macro, pooled, rt), encoding="utf-8")

    git = git_state()["commit"]
    rti = rt.set_index(["crop_id", "method"])
    rows = []
    for r in per[per.method.str.startswith("cellpose:")].itertuples():
        tag = r.method.split(":", 1)[1]
        t = rti.loc[(r.crop_id, r.method)] if (r.crop_id, r.method) in rti.index else None
        rj = RUNS / tag / r.crop_id / ("run.json" if r.status == "ok" else "FAILED.json")
        rec = json.loads(rj.read_text(encoding="utf-8"))
        rows.append({"experiment_id": f"cellpose_{tag}/{r.crop_id}", "method": "cellpose", "model": rec["model"],
                     "volume": r.scan, "crop": r.crop_id,
                     "device": t.device if t is not None else rec["environment"]["device"],
                     "parameters": f"size_filter={r.size_filter};eval=evaluate_v2:{MODE};split={r.split};"
                                   f"tag={tag};overrides={json.dumps(rec.get('overrides', {}), sort_keys=True)};"
                                   f"peak_vram_reserved_mb={None if t is None else t.peak_vram_reserved_mb}",
                     "runtime_seconds": None if t is None else t.seconds,
                     "peak_memory_mb": None if t is None else t.peak_host_mb,
                     "tp": r.tp_iou050, "fp": r.fp_iou050, "fn": r.fn_iou050, "precision": r.precision_iou050,
                     "recall": r.recall_iou050, "f1_025": r.f1_iou025, "f1_050": r.f1_iou050, "f1_075": r.f1_iou075,
                     "ap50": r.ap_iou050, "ap75": r.ap_iou075, "mean_matched_iou": r.mean_matched_iou050,
                     "n_pred": r.prediction_count_evaluable, "n_gt": r.valid_gt_count, "git_commit": git,
                     "run_json": str(rj)})
    upsert(REPO / "results" / "all_runs.csv", rows, ["experiment_id", "method", "parameters"])
    print((out / "tables.md").read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
