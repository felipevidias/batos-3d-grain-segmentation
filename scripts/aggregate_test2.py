#!/usr/bin/env python3
"""Phase 5f step 4: frozen / fixflat / v3 / v3.1 side by side on held-out TEST SET 2, plus the test-set-1 v3 numbers
(results/test_v3, unchanged) for reference.

Inputs : C:\\BA_TOS_DATA\\runs\\test2\\eval\\<crop>.{csv,json}      (frozen, fixflat; scripts/run_test2.py)
         C:\\BA_TOS_DATA\\runs\\test2\\eval_v3\\<crop>.{csv,json}   (v3)
         C:\\BA_TOS_DATA\\runs\\test2\\eval_v31\\<crop>.{csv,json}  (v3.1)
Outputs: results/test2/ per_crop_all_modes.csv, crop_info.csv, crop_macro.csv, pooled.csv, per_crop_primary.csv,
         scale.csv, runtime.csv, forensics_stages.csv, flags.csv, tables.md; rows in results/all_runs.csv
         (experiment_id test2_<method>/<crop>).
Scopes, metrics and the pre-registered first-stage rule are those of scripts/aggregate_test_v3.py (imported). Added:
seed density (MIN candidates and ToS markers per GT grain incl. face-touching; GT counts = diagnosis only).
Flags: v3.1 F1@0.50 < 0.90 (ToS or sel, markers/dominant, no filter); v3 flags listed as well.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
import aggregate_test_v3 as agg  # noqa: E402
from batos_v2.extinction_run import upsert  # noqa: E402
from batos_v2.run_crop import git_state  # noqa: E402

ROOT = Path(r"C:\BA_TOS_DATA\runs\test2")
OUT = REPO / "results" / "test2"
T1 = REPO / "results" / "test_v3"
METHODS = ["frozen", "fixflat", "v3", "v3.1"]
ORDER = {m: i for i, m in enumerate(METHODS)}
N_CROPS = 52
RUN_DIR = {"v3": "v3", "v3.1": "v31"}
FLAG_F1 = agg.FLAG_F1


def read(sub: str, ext: str) -> list:
    return sorted((ROOT / sub).glob(f"*.{ext}"))


def load():
    parts = {s: pd.concat([pd.read_csv(p, dtype={"size_filter": str}) for p in read(s, "csv")])
             for s in ("eval", "eval_v3", "eval_v31")}
    ids = {s: sorted(d.crop_id.unique()) for s, d in parts.items()}
    if len({tuple(v) for v in ids.values()}) != 1 or len(ids["eval"]) != N_CROPS:
        raise RuntimeError(f"crop sets differ or incomplete: { {s: len(v) for s, v in ids.items()} }")
    if set(parts["eval_v31"].method) != {"v3.1"} or set(parts["eval_v3"].method) != {"v3"}:
        raise RuntimeError("unexpected method labels in v3 / v3.1 evaluations")
    per = pd.concat([parts["eval"]] + [parts[s].drop(columns=["size_filter_voxels"]) for s in ("eval_v3", "eval_v31")],
                    ignore_index=True)
    js = {s: pd.DataFrame([json.loads(p.read_text(encoding="utf-8")) for p in read(s, "json")])
          for s in ("eval", "eval_v3", "eval_v31")}
    for k in ("frozen_failed", "fixflat_failed"):
        if k not in js["eval"]:
            js["eval"][k] = None
    drop = ["specimen", "scan", "test_class", "eval_seconds"]
    i3 = js["eval_v3"].drop(columns=drop)
    i31 = js["eval_v31"].drop(columns=drop).rename(columns=lambda c: c.replace("v3_", "v31_", 1))
    info = js["eval"].merge(i3, on="crop_id").merge(i31, on="crop_id")
    return per, info


class _Root:
    """Stand-in for agg.ROOT so that agg.v3_stages(crop) reads <ROOT>/<run dir of the method>/<crop>."""

    def __init__(self, sub: str):
        self.sub = sub

    def __truediv__(self, name):
        return ROOT / (self.sub if name == "v3" else name)


def stages(method: str, crop: str, x: pd.Series) -> dict:
    """GT-free stage quantities of one v3 / v3.1 run (agg.v3_stages, unchanged, on the method's run directory)."""
    old = agg.ROOT
    agg.ROOT = _Root(RUN_DIR[method])
    try:
        r = agg.v3_stages(crop)
    finally:
        agg.ROOT = old
    gt_all = int(x.valid_gt_count + x.n_gt_excluded)
    r |= {"gt_valid": int(x.valid_gt_count), "gt_all": gt_all,
          "candidates_per_grain": (r.get("min_candidates") or np.nan) / gt_all,
          "markers_per_grain": (r.get("tos_markers") or np.nan) / gt_all}
    return r


def md_macro(t, L, title):
    L.append(f"\n### {title}\n")
    L.append("| scope | route | method | n | F1@0.50 mean ± SD | median | min (crop) | # < 0.90 | F1@0.75 | F1@0.25 |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in t.itertuples():
        L.append(f"| {r.scope} | {r.route} | {r.method} | {r.n_crops} | {r.f1_050_mean:.3f} ± {r.f1_050_sd:.3f} "
                 f"| {r.f1_050_median:.3f} | {r.f1_050_min:.3f} ({r.f1_050_argmin}) | {r.n_below_090_050} "
                 f"| {r.f1_075_mean:.3f} | {r.f1_025_mean:.3f} |")


def md_pooled(t, L, title):
    L.append(f"\n### {title}\n")
    L.append("| scope | route | method | GT | pred | TP | FP | FN | P@0.50 | R@0.50 | F1@0.50 | F1@0.75 |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in t.itertuples():
        L.append(f"| {r.scope} | {r.route} | {r.method} | {r.valid_gt} | {r.pred_eval} | {r.tp_050} | {r.fp_050} "
                 f"| {r.fn_050} | {r.precision_050:.3f} | {r.recall_050:.3f} | {r.f1_050:.3f} | {r.f1_075:.3f} |")


def sort(t):
    t = t.copy()
    t["o"] = t.method.map(ORDER)
    return t.sort_values(["scope", "route", "o"], kind="stable")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    per, info = load()
    per.to_csv(OUT / "per_crop_all_modes.csv", index=False)
    info.to_csv(OUT / "crop_info.csv", index=False)
    macro, pooled = agg.aggregate(per)
    macro.to_csv(OUT / "crop_macro.csv", index=False)
    pooled.to_csv(OUT / "pooled.csv", index=False)

    prim = per[(per["mode"] == "markers_dominant") & (per.size_filter == "none")]
    wide = prim.pivot_table(index="crop_id", columns=["method", "route"], values="f1_iou050")
    wide.columns = [f"{a.replace('.', '')}_{b[7:]}_f1_050" for a, b in wide.columns]
    wide = wide.join(info.set_index("crop_id"))
    wide.to_csv(OUT / "per_crop_primary.csv")

    sc = info[["crop_id", "specimen", "test_class", "gt_inscribed_r_median", "r_est_e13", "v3_r_granulo",
               "v31_r_granulo", "v3_openrec_radius"]].copy()
    sc["ratio_granulo"] = sc.v3_r_granulo / sc.gt_inscribed_r_median
    sc["ratio_e13"] = sc.r_est_e13 / sc.gt_inscribed_r_median
    sc.to_csv(OUT / "scale.csv", index=False)
    ok = sc.dropna(subset=["v3_r_granulo", "gt_inscribed_r_median"])
    ok13 = sc.dropna(subset=["r_est_e13", "gt_inscribed_r_median"])
    corr = float(np.corrcoef(ok.v3_r_granulo, ok.gt_inscribed_r_median)[0, 1])
    corr13 = float(np.corrcoef(ok13.r_est_e13, ok13.gt_inscribed_r_median)[0, 1])

    rt = pd.DataFrame({"crop_id": info.crop_id, "specimen": info.specimen, "frozen_s": info.frozen_seconds,
                       "fixflat_s": info.fixflat_seconds_end_to_end, "v3_s": info.v3_seconds, "v31_s": info.v31_seconds,
                       "v3_peak_mb": info.v3_peak_mb, "v31_peak_mb": info.v31_peak_mb,
                       "openrec_impl": info.v3_openrec_impl})
    rt.to_csv(OUT / "runtime.csv", index=False)

    frows = []
    ii = info.set_index("crop_id")
    for crop in info.crop_id:
        i = ii.loc[crop]
        for method in ("v3", "v3.1"):
            p = "v3" if method == "v3" else "v31"
            x = prim[(prim.method == method) & (prim.crop_id == crop) & (prim.route == "BA-TOS_sel")].iloc[0]
            t = prim[(prim.method == method) & (prim.crop_id == crop) & (prim.route == "BA-TOS_ToS")].iloc[0]
            r = {"crop_id": crop, "method": method, "specimen": i.specimen, "test_class": i.test_class,
                 **stages(method, crop, x), "r_granulo": i[f"{p}_r_granulo"], "gt_r_in": i.gt_inscribed_r_median,
                 "openrec_radius": i[f"{p}_openrec_radius"], "stage03_thr": i[f"{p}_stage03_thr"],
                 "gate": i[f"{p}_gate_triggered"], "mv": i[f"{p}_mv"], "qv": i[f"{p}_qv"],
                 "valid_gt": int(x.valid_gt_count), "precision": float(x.precision_iou050),
                 "recall": float(x.recall_iou050), "f1_tos": float(t.f1_iou050), "f1_sel": float(x.f1_iou050),
                 "fixflat_f1_sel": wide.loc[crop, "fixflat_sel_f1_050"], "frozen_f1_sel": wide.loc[crop, "frozen_sel_f1_050"]}
            r["r_ratio"] = r["r_granulo"] / r["gt_r_in"] if r["r_granulo"] and r["gt_r_in"] else None
            r["flag"] = bool(min(r["f1_tos"], r["f1_sel"]) < FLAG_F1)
            r["first_stage"] = agg.first_stage(r) if r["flag"] else ""
            frows.append(r)
    fr = pd.DataFrame(frows)
    fr.to_csv(OUT / "forensics_stages.csv", index=False)
    flags = fr[fr.flag]
    flags.to_csv(OUT / "flags.csv", index=False)

    L = ["# Test set 2 tables (generated by scripts/aggregate_test2.py)"]
    for flt in ("none", "500", "own"):
        note = " (v3 / v3.1 only: 500 s^3)" if flt == "own" else ""
        for kind, tab, fn in (("Crop-macro F1", macro, md_macro), ("Pooled objects", pooled, md_pooled)):
            t = tab[(tab["mode"] == "markers_dominant") & (tab.size_filter == flt)]
            if not t.empty:
                fn(sort(t), L, f"{kind} (markers/dominant = seed-ignore), size filter {flt}{note}")
    t = macro[(macro["mode"] == "no_markers_fraction050") & macro.scope.isin(["all", "EFR_internal", "external"])]
    L.append("\n### No-marker mode (own-voxel support, overlap-ignore 0.50), crop-macro F1@0.50\n")
    L.append("| scope | route | filter | method | mean ± SD | min |")
    L.append("|---|---|---|---|---|---|")
    for r in sort(t).sort_values(["scope", "route", "size_filter", "o"]).itertuples():
        L.append(f"| {r.scope} | {r.route} | {r.size_filter} | {r.method} | {r.f1_050_mean:.3f} ± {r.f1_050_sd:.3f} "
                 f"| {r.f1_050_min:.3f} |")

    L.append("\n### Test set 1 vs test set 2 (crop-macro and pooled F1@0.50, BA-TOS_sel, seed-ignore, no filter)\n")
    m1 = pd.read_csv(T1 / "crop_macro.csv", dtype={"size_filter": str})
    p1 = pd.read_csv(T1 / "pooled.csv", dtype={"size_filter": str})
    sel = lambda d: d[(d["mode"] == "markers_dominant") & (d.size_filter == "none") & (d.route == "BA-TOS_sel")]  # noqa: E731
    L.append("| scope | set | frozen | fixflat | v3 | v3.1 | pooled frozen | pooled fixflat | pooled v3 | pooled v3.1 |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for scope in ["all", "EFR_internal", "external", "EFR", "EFS", "EMR", "EMS", "EOR", "EOS"]:
        for name, m, p in (("test 1", sel(m1), sel(p1)), ("test 2", sel(macro), sel(pooled))):
            mm = m[m.scope == scope].set_index("method")
            pp = p[p.scope == scope].set_index("method")
            cell = lambda d, k, c: f"{d.loc[k, c]:.3f}" if k in d.index else "—"  # noqa: E731
            L.append(f"| {scope} | {name} | " + " | ".join(cell(mm, k, "f1_050_mean") for k in METHODS) + " | "
                     + " | ".join(cell(pp, k, "f1_050") for k in METHODS) + " |")

    L.append("\n### Gate and failures\n")
    L.append("| specimen | n | gate fires frozen / fixflat / v3 / v3.1 | crash frozen / fixflat | v3 crash / empty | v3.1 crash / empty |")
    L.append("|---|---|---|---|---|---|")
    for spec, g in list(info.groupby("specimen")) + [("all", info)]:
        fire = lambda c: int((g[c].astype(str) == "True").sum())  # noqa: E731
        L.append(f"| {spec} | {len(g)} | {fire('frozen_gate_triggered')} / {fire('fixflat_gate_triggered')} / "
                 f"{fire('v3_gate_triggered')} / {fire('v31_gate_triggered')} | {int(g.frozen_failed.notna().sum())} / "
                 f"{int(g.fixflat_failed.notna().sum())} | {int(g.v3_failed.notna().sum())} / {int(g.v3_empty.notna().sum())} "
                 f"| {int(g.v31_failed.notna().sum())} / {int(g.v31_empty.notna().sum())} |")

    L.append(f"\n### Scale: r_granulo vs GT median inscribed radius (GT diagnostic only); Pearson r = {corr:.3f} "
             f"(e13: {corr13:.3f}, n = {len(ok13)}), n = {len(ok)}\n")
    L.append("| specimen | n | GT r_in | r_granulo median [min–max] | ratio | r_e13 median | ratio |")
    L.append("|---|---|---|---|---|---|---|")
    for spec, g in sc.groupby("specimen"):
        L.append(f"| {spec} | {len(g)} | {g.gt_inscribed_r_median.median():.2f} | {g.v3_r_granulo.median():g} "
                 f"[{g.v3_r_granulo.min():g}–{g.v3_r_granulo.max():g}] | {g.ratio_granulo.median():.2f} "
                 f"| {g.r_est_e13.median():.2f} | {g.ratio_e13.median():.2f} |")

    L.append("\n### Seed density (v3 vs v3.1; per GT grain incl. face-touching, GT counts = diagnosis only)\n")
    L.append("| specimen | class | n | v3 cand./grain | v3.1 cand./grain | v3 markers/grain | v3.1 markers/grain "
             "| v3 F1 median | v3.1 F1 median |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    w = fr.pivot_table(index=["crop_id", "specimen", "test_class"], columns="method",
                       values=["candidates_per_grain", "markers_per_grain", "f1_sel"]).reset_index()
    for (spec, cls), g in w.groupby(["specimen", "test_class"]):
        L.append(f"| {spec} | {cls} | {len(g)} | {g[('candidates_per_grain', 'v3')].median():.2f} "
                 f"| {g[('candidates_per_grain', 'v3.1')].median():.2f} | {g[('markers_per_grain', 'v3')].median():.2f} "
                 f"| {g[('markers_per_grain', 'v3.1')].median():.2f} | {g[('f1_sel', 'v3')].median():.3f} "
                 f"| {g[('f1_sel', 'v3.1')].median():.3f} |")

    L.append("\n### Runtime per crop (s, end to end; not contention-controlled)\n")
    L.append("| specimen | n | frozen | fixflat | v3 | v3.1 | v3.1 max | v3.1 peak MB max |")
    L.append("|---|---|---|---|---|---|---|---|")
    for spec, g in list(rt.groupby("specimen")) + [("all", rt)]:
        L.append(f"| {spec} | {len(g)} | {g.frozen_s.median():.0f} | {g.fixflat_s.median():.0f} | {g.v3_s.median():.0f} "
                 f"| {g.v31_s.median():.0f} | {g.v31_s.max():.0f} | {g.v31_peak_mb.max():.0f} |")
    L.append(f"\nTotals (h): frozen {rt.frozen_s.sum() / 3600:.1f}, fixflat {rt.fixflat_s.sum() / 3600:.1f}, "
             f"v3 {rt.v3_s.sum() / 3600:.1f}, v3.1 {rt.v31_s.sum() / 3600:.1f}\n")

    fmt = lambda v, f="{:g}": "—" if v is None or (isinstance(v, float) and np.isnan(v)) else f.format(v)  # noqa: E731
    for method in ("v3.1", "v3"):
        fl = flags[flags.method == method]
        L.append(f"\n### Flagged crops: {method} F1@0.50 < {FLAG_F1} (ToS or sel) — {len(fl)} crops, first-stage forensics\n")
        L.append("| crop | class | ToS / sel | fixflat sel | frozen sel | r / GT r_in | openrec r | thr03 "
                 "| seeds kept / total | MIN cand. | cand. / grain | ToS markers / grain | gate (mv/qv) | P / R | first stage (frozen rule) |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in fl.to_dict("records"):
            L.append(f"| {r['crop_id']} | {r['test_class']} | {r['f1_tos']:.3f} / {r['f1_sel']:.3f} | {r['fixflat_f1_sel']:.3f} "
                     f"| {r['frozen_f1_sel']:.3f} | {fmt(r['r_granulo'])} / {fmt(r['gt_r_in'], '{:.1f}')} "
                     f"| {fmt(r['openrec_radius'])} | {fmt(r['stage03_thr'])} | {fmt(r.get('bg_keep'))} / "
                     f"{fmt(r.get('seed_components'))} | {fmt(r.get('min_candidates'))} "
                     f"| {fmt(r['candidates_per_grain'], '{:.2f}')} | {fmt(r['markers_per_grain'], '{:.2f}')} "
                     f"| {r['gate']} ({fmt(r['mv'])}/{fmt(r['qv'])}) | {r['precision']:.2f} / {r['recall']:.2f} "
                     f"| {r['first_stage']} |")

    L.append("\n### Per crop (markers/dominant, no filter, F1@0.50)\n")
    L.append("| crop | class | frozen ToS / sel | fixflat ToS / sel | v3 ToS / sel | v3.1 ToS / sel | r_granulo / GT r_in |")
    L.append("|---|---|---|---|---|---|---|")
    for cid, r in wide.iterrows():
        L.append(f"| {cid} | {r.test_class} | {r.frozen_ToS_f1_050:.3f} / {r.frozen_sel_f1_050:.3f} | "
                 f"{r.fixflat_ToS_f1_050:.3f} / {r.fixflat_sel_f1_050:.3f} | {r.v3_ToS_f1_050:.3f} / {r.v3_sel_f1_050:.3f} | "
                 f"{r.v31_ToS_f1_050:.3f} / {r.v31_sel_f1_050:.3f} | {fmt(r.v3_r_granulo)} / {fmt(r.gt_inscribed_r_median, '{:.1f}')} |")
    (OUT / "tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    git = git_state()["commit"]
    model = {"frozen": "batos-paper-v2-frozen", "fixflat": "batos-paper-v2-fixflat", "v3": "batos-paper-v3-scalefree",
             "v3.1": "batos-paper-v3.1"}
    rows = []
    for r in per.itertuples():
        i = ii.loc[r.crop_id]
        sec = {"frozen": i.frozen_seconds, "fixflat": i.fixflat_seconds_end_to_end, "v3": i.v3_seconds,
               "v3.1": i.v31_seconds}[r.method]
        rows.append({"experiment_id": f"test2_{r.method}/{r.crop_id}", "method": r.route, "model": model[r.method],
                     "volume": r.scan, "crop": r.crop_id, "device": "cpu (AMD Ryzen 7 5700X)",
                     "parameters": f"size_filter={r.size_filter};eval=evaluate_v2:{r.mode};split=test2",
                     "runtime_seconds": sec, "tp": r.tp_iou050, "fp": r.fp_iou050, "fn": r.fn_iou050,
                     "precision": r.precision_iou050, "recall": r.recall_iou050, "f1_025": r.f1_iou025,
                     "f1_050": r.f1_iou050, "f1_075": r.f1_iou075, "ap50": r.ap_iou050, "ap75": r.ap_iou075,
                     "mean_matched_iou": r.mean_matched_iou050, "n_pred": r.prediction_count_evaluable,
                     "n_gt": r.valid_gt_count, "git_commit": git, "run_json": str(ROOT)})
    upsert(REPO / "results" / "all_runs.csv", rows, ["experiment_id", "method", "parameters"])
    print((OUT / "tables.md").read_text(encoding="utf-8"))
    for m in ("v3.1", "v3"):
        print(f"flags {m} (F1@0.50 < {FLAG_F1}): {list(flags[flags.method == m].crop_id)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
