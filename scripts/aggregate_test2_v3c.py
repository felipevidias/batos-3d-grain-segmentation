#!/usr/bin/env python3
"""Phase 5i: v3 vs v3c on held-out TEST SET 2 (scripts/run_test2_v3c.py), per specimen (EOR / EOS emphasised).

Inputs : C:\\BA_TOS_DATA\\runs\\test2_v3c\\eval\\<crop>.{csv,json}; data/metadata/crops_test2_v3c.json (native sizes).
Outputs: results/test2_v3c/{per_crop_all.csv, crop_info.csv, crop_macro.csv, pooled.csv, per_crop_primary.csv,
         native_sizes.csv, seed_boundary.csv, tables.md}; rows in results/all_runs.csv (experiment_id
         test2_v3c_<method>/<crop>).
Primary: BA-TOS_sel, markers_dominant (seed-ignore), no size filter, F1@0.50; v3 with valid GT 'fixed1000' (its own
protocol on the 200^3 probe), v3c with 'canon1000' (>= round(1000 (r/10)^3) native voxels on its window). The other rule
is reported for both. Crop-macro (mean of per-crop F1) and pooled (sum of TP/FP/FN) are separate tables.
Scopes: each specimen, EFR_internal, external (5 specimens), all 52.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.extinction_run import upsert  # noqa: E402
from batos_v2.run_crop import git_state  # noqa: E402

ROOT = Path(r"C:\BA_TOS_DATA\runs\test2_v3c")
OUT = REPO / "results" / "test2_v3c"
N_CROPS = 52
SPECIMENS = ["EFR", "EFS", "EMR", "EMS", "EOR", "EOS"]
KEYS = ["method", "gt_rule", "route", "size_filter", "mode"]
PRIMARY = {"v3": "fixed1000", "v3c": "canon1000"}


def scopes(df: pd.DataFrame):
    for s in SPECIMENS:
        yield s, df[df.specimen == s]
    yield "EFR_internal", df[df.specimen == "EFR"]
    yield "external", df[df.specimen != "EFR"]
    yield "all", df


def macro(per: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scope, d in scopes(per):
        for k, g in d.groupby(KEYS, dropna=False):
            row = dict(zip(KEYS, k)) | {"scope": scope, "n_crops": g.crop_id.nunique()}
            for t in ("025", "050", "075"):
                x = g[f"f1_iou{t}"]
                row |= {f"f1_{t}_mean": x.mean(), f"f1_{t}_sd": x.std(ddof=1), f"f1_{t}_median": x.median(),
                        f"f1_{t}_min": x.min()}
            rows.append(row)
    return pd.DataFrame(rows)


def pooled(per: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scope, d in scopes(per):
        for k, g in d.groupby(KEYS, dropna=False):
            row = dict(zip(KEYS, k)) | {"scope": scope, "n_crops": g.crop_id.nunique()}
            for t in ("025", "050", "075"):
                tp, fp, fn = (int(g[f"{m}_iou{t}"].sum()) for m in ("tp", "fp", "fn"))
                p = tp / (tp + fp) if tp + fp else 0.0
                r = tp / (tp + fn) if tp + fn else 0.0
                row |= {f"tp_{t}": tp, f"fp_{t}": fp, f"fn_{t}": fn, f"precision_{t}": p, f"recall_{t}": r,
                        f"f1_{t}": 2 * p * r / (p + r) if p + r else 0.0}
            rows.append(row)
    return pd.DataFrame(rows)


def fmt(x, f="{:.3f}"):
    return "–" if x is None or (isinstance(x, float) and np.isnan(x)) else f.format(x)


def main() -> int:
    csvs = sorted((ROOT / "eval").glob("*.csv"))
    per = pd.concat([pd.read_csv(p, dtype={"size_filter": str}) for p in csvs], ignore_index=True)
    info = pd.DataFrame([json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "eval").glob("*.json"))])
    if per.crop_id.nunique() != N_CROPS or len(info) != N_CROPS:
        raise RuntimeError(f"incomplete: {per.crop_id.nunique()} crops evaluated, {len(info)} info files (need {N_CROPS})")
    win = pd.DataFrame(json.loads((REPO / "data" / "metadata" / "crops_test2_v3c.json").read_text(encoding="utf-8"))["crops"])
    OUT.mkdir(parents=True, exist_ok=True)
    per.to_csv(OUT / "per_crop_all.csv", index=False)
    info.to_csv(OUT / "crop_info.csv", index=False)
    cm, po = macro(per), pooled(per)
    cm.to_csv(OUT / "crop_macro.csv", index=False)
    po.to_csv(OUT / "pooled.csv", index=False)
    sizes = win[["crop_id", "specimen", "test_class", "r_probe", "native_side", "origin_zyx", "shift_zyx",
                 "centred_exactly", "min_gt_voxels", "valid_gt"]]
    sizes.to_csv(OUT / "native_sizes.csv", index=False)
    prim = per[(per.route == "BA-TOS_sel") & (per["mode"] == "markers_dominant") & (per.size_filter == "none")]
    wide = prim.pivot_table(index=["crop_id", "specimen", "test_class"], columns=["method", "gt_rule"],
                            values="f1_iou050").reset_index()
    wide.columns = ["_".join(c for c in col if c) if isinstance(col, tuple) else col for col in wide.columns]
    wide = wide.merge(info[["crop_id", "r_probe", "native_side", "v3_r_granulo", "v3_gate", "v3c_gate", "v3_failed",
                            "v3_empty", "v3c_failed", "bnd_n_common", "bnd_v3c_miou_common", "bnd_v3_miou_common",
                            "bnd_v3c_volerr_median", "bnd_v3_volerr_median", "bnd_v3c_abs_volerr_median",
                            "bnd_v3_abs_volerr_median", "v3_seed_markers_per_valid_grain",
                            "v3c_seed_markers_per_valid_grain", "v3_seed_valid_split", "v3c_seed_valid_split"]],
                      on="crop_id")
    wide.to_csv(OUT / "per_crop_primary.csv", index=False)
    sb_cols = [c for c in info.columns if c.startswith(("bnd_", "v3_seed_", "v3c_seed_"))]
    sb = []
    for scope, d in scopes(info):
        row = {"scope": scope, "n_crops": len(d)}
        for c in sb_cols:
            x = pd.to_numeric(d[c], errors="coerce")
            row[c + ("_sum" if c.endswith(("n_common", "_split", "_unseeded", "matched_v3", "matched_v3c")) else "_median")] = \
                x.sum() if c.endswith(("n_common", "_split", "_unseeded", "matched_v3", "matched_v3c")) else x.median()
        sb.append(row)
    sb = pd.DataFrame(sb)
    sb.to_csv(OUT / "seed_boundary.csv", index=False)

    def cmv(method, scope, rule=None, route="BA-TOS_sel", flt="none", mode="markers_dominant", col="f1_050_mean"):
        r = cm[(cm.method == method) & (cm.gt_rule == (rule or PRIMARY[method])) & (cm.route == route)
               & (cm.size_filter == flt) & (cm["mode"] == mode) & (cm.scope == scope)]
        return float(r[col].iloc[0]) if len(r) else float("nan")

    def pov(method, scope, rule=None, t="050"):
        r = po[(po.method == method) & (po.gt_rule == (rule or PRIMARY[method])) & (po.route == "BA-TOS_sel")
               & (po.size_filter == "none") & (po["mode"] == "markers_dominant") & (po.scope == scope)]
        return float(r[f"f1_{t}"].iloc[0]) if len(r) else float("nan")

    sc = SPECIMENS + ["EFR_internal", "external", "all"]
    L = ["# Test set 2 — v3 vs v3c (Phase 5i, run once after the v3c freeze)", "",
         "Primary: BA-TOS_sel, seed-ignore (markers_dominant), no size filter. v3 on the 200³ probe with valid GT ≥ 1000 "
         "voxels; v3c on its canonical-FOV window (native side round(200 r/10)) with valid GT ≥ round(1000 (r/10)³) "
         "native voxels. Each method against its own window's GT. The cross rule is shown in parentheses.", "",
         "## Crop-macro F1@0.50 (mean ± SD, min) — EOR / EOS emphasised", "",
         "| scope | n | v3 | v3c | Δ (v3c − v3) | v3 (canon rule) | v3c (fixed rule) |", "|---|---|---|---|---|---|---|"]
    for s in sc:
        n = int(cm[(cm.scope == s)].n_crops.max()) if (cm.scope == s).any() else 0
        a, b = cmv("v3", s), cmv("v3c", s)
        bold = "**" if s in ("EOR", "EOS") else ""
        L.append(f"| {bold}{s}{bold} | {n} | {fmt(a)} ± {fmt(cmv('v3', s, col='f1_050_sd'))} (min {fmt(cmv('v3', s, col='f1_050_min'))}) "
                 f"| {fmt(b)} ± {fmt(cmv('v3c', s, col='f1_050_sd'))} (min {fmt(cmv('v3c', s, col='f1_050_min'))}) "
                 f"| {bold}{b - a:+.3f}{bold} | {fmt(cmv('v3', s, 'canon1000'))} | {fmt(cmv('v3c', s, 'fixed1000'))} |")
    L += ["", "## Crop-macro F1@0.25 / F1@0.75 and ToS route (primary rules)", "",
          "| scope | v3 sel .25 / .75 | v3c sel .25 / .75 | v3 ToS .50 | v3c ToS .50 | v3 sel own-filter .50 | v3c sel own-filter .50 |",
          "|---|---|---|---|---|---|---|"]
    for s in sc:
        L.append(f"| {s} | {fmt(cmv('v3', s, col='f1_025_mean'))} / {fmt(cmv('v3', s, col='f1_075_mean'))} "
                 f"| {fmt(cmv('v3c', s, col='f1_025_mean'))} / {fmt(cmv('v3c', s, col='f1_075_mean'))} "
                 f"| {fmt(cmv('v3', s, route='BA-TOS_ToS'))} | {fmt(cmv('v3c', s, route='BA-TOS_ToS'))} "
                 f"| {fmt(cmv('v3', s, flt='own'))} | {fmt(cmv('v3c', s, flt='own'))} |")
    L += ["", "## Pooled F1 (sum of TP/FP/FN over crops; separate from crop-macro)", "",
          "| scope | v3 .25 / .50 / .75 | v3c .25 / .50 / .75 |", "|---|---|---|"]
    for s in sc:
        L.append(f"| {s} | {' / '.join(fmt(pov('v3', s, t=t)) for t in ('025', '050', '075'))} "
                 f"| {' / '.join(fmt(pov('v3c', s, t=t)) for t in ('025', '050', '075'))} |")
    L += ["", "## Boundary cost of the resampling (sel; grains valid in both windows and matched IoU > 0.5 by both)", "",
          "Per-crop mean matched IoU and median relative volume error; the table gives the median over crops and the "
          "number of common grains.", "",
          "| scope | common grains | mIoU v3 | mIoU v3c | median vol. err. v3 | v3c | median |vol. err.| v3 | v3c |",
          "|---|---|---|---|---|---|---|---|"]
    for r in sb.itertuples():
        L.append(f"| {r.scope} | {int(r.bnd_n_common_sum)} | {fmt(r.bnd_v3_miou_common_median)} | {fmt(r.bnd_v3c_miou_common_median)} "
                 f"| {fmt(r.bnd_v3_volerr_median_median, '{:+.3f}')} | {fmt(r.bnd_v3c_volerr_median_median, '{:+.3f}')} "
                 f"| {fmt(r.bnd_v3_abs_volerr_median_median)} | {fmt(r.bnd_v3c_abs_volerr_median_median)} |")
    L += ["", "## Seeds (sel markers per valid grain, canon rule; GT = diagnosis only)", "",
          "| scope | v3 markers/grain (median over crops) | v3c | v3 split grains | v3c split grains |", "|---|---|---|---|---|"]
    for r in sb.itertuples():
        L.append(f"| {r.scope} | {fmt(r.v3_seed_markers_per_valid_grain_median)} | {fmt(r.v3c_seed_markers_per_valid_grain_median)} "
                 f"| {int(r.v3_seed_valid_split_sum)} | {int(r.v3c_seed_valid_split_sum)} |")
    L += ["", "## Native window sizes used by v3c", "", "| specimen | r_probe values | native sides | shifted windows |",
          "|---|---|---|---|"]
    for s in SPECIMENS:
        d = sizes[sizes.specimen == s]
        L.append(f"| {s} | {sorted(set(d.r_probe))} | {sorted(set(d.native_side))} | "
                 f"{int((~d.centred_exactly.astype(bool)).sum())} / {len(d)} |")
    L += ["", "## Per crop (sel, seed-ignore, no filter, F1@0.50)", "",
          "| crop | class | r | S | v3 | v3c | Δ | gate v3 / v3c | mIoU v3 / v3c |", "|---|---|---|---|---|---|---|---|---|"]
    for r in wide.sort_values(["specimen", "crop_id"]).itertuples():
        st3 = "CRASH" if isinstance(r.v3_failed, str) else ("EMPTY" if isinstance(r.v3_empty, str) or r.v3_empty is True else "")
        stc = "CRASH" if isinstance(r.v3c_failed, str) else ""
        L.append(f"| {r.crop_id} | {r.test_class} | {r.r_probe:g} | {r.native_side} | {r.v3_fixed1000:.3f} {st3} "
                 f"| {r.v3c_canon1000:.3f} {stc} | {r.v3c_canon1000 - r.v3_fixed1000:+.3f} | {r.v3_gate} / {r.v3c_gate} "
                 f"| {fmt(r.bnd_v3_miou_common)} / {fmt(r.bnd_v3c_miou_common)} |")
    (OUT / "tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    git = git_state()["commit"]
    ii, rows = info.set_index("crop_id"), []
    for r in per.itertuples():
        i = ii.loc[r.crop_id]
        rows.append({"experiment_id": f"test2_v3c_{r.method}/{r.crop_id}", "method": r.route,
                     "model": "batos-paper-v3c" if r.method == "v3c" else "batos-paper-v3-scalefree",
                     "volume": r.scan, "crop": r.crop_id, "device": "cpu (AMD Ryzen 7 5700X)",
                     "parameters": f"size_filter={r.size_filter};gt_rule={r.gt_rule};eval=eval_lean:{r.mode};split=test2;"
                                   f"r_granulo={i.r_probe};native_side={i.native_side if r.method == 'v3c' else 200}",
                     "runtime_seconds": i.v3c_seconds if r.method == "v3c" else None,
                     "peak_memory_mb": i.v3c_peak_mb if r.method == "v3c" else None, "tp": r.tp_iou050,
                     "fp": r.fp_iou050, "fn": r.fn_iou050, "precision": r.precision_iou050, "recall": r.recall_iou050,
                     "f1_025": r.f1_iou025, "f1_050": r.f1_iou050, "f1_075": r.f1_iou075, "ap50": r.ap_iou050,
                     "ap75": r.ap_iou075, "mean_matched_iou": r.mean_matched_iou050,
                     "n_pred": r.prediction_count_evaluable, "n_gt": r.valid_gt_count, "git_commit": git,
                     "run_json": str(ROOT / "eval" / f"{r.crop_id}.json")})
    upsert(REPO / "results" / "all_runs.csv", rows, ["experiment_id", "method", "parameters"])
    print((OUT / "tables.md").read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
