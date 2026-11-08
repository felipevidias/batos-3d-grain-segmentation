#!/usr/bin/env python3
"""Phase 5e step 4: frozen / fixflat / v3 side by side on the 52 test crops.

Inputs : C:\\BA_TOS_DATA\\runs\\test_protocol\\eval\\<crop>.{csv,json}    (Phase 1f: frozen, fixflat; reused, not re-run)
         C:\\BA_TOS_DATA\\runs\\test_protocol\\eval_v3\\<crop>.{csv,json} (scripts/run_test_v3.py)
Outputs: results/test_v3/ per_crop_all_modes.csv, crop_info.csv, crop_macro.csv, pooled.csv, per_crop_primary.csv,
         scale.csv, runtime.csv, forensics_stages.csv, flags.csv, tables.md; rows in results/all_runs.csv
         (experiment_id test_v3_v3/<crop>; frozen/fixflat rows already exist as test_protocol_<method>/<crop>).
Scopes : each specimen, EFR_internal, external (5 specimens), all 52. Crop-macro and pooled kept separate.
Size filter: none (primary) and 500 for all methods; v3 also 'own' (500 s^3).
Flags  : v3 F1@0.50 < 0.90 (ToS or sel, markers/dominant, no filter) with first-stage forensics (GT-free stage
         quantities from the kept intermediates; GT only via F1, valid-GT count and inscribed radius, diagnosis only).
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.extinction_run import upsert  # noqa: E402
from batos_v2.run_crop import git_state  # noqa: E402

ROOT = Path(r"C:\BA_TOS_DATA\runs\test_protocol")
OUT = REPO / "results" / "test_v3"
KEYS = ["method", "route", "size_filter", "mode"]
METHODS = ["frozen", "fixflat", "v3"]
FLAG_F1 = 0.90
SCALE_TOL = 0.30  # forensics rule: |r / GT r_in - 1| > 0.30 marks the scale stage (diagnostic only)


def scopes(df):
    yield "all", df
    yield "EFR_internal", df[df.specimen == "EFR"]
    yield "external", df[df.specimen != "EFR"]
    for s in sorted(df.specimen.unique()):
        yield s, df[df.specimen == s]


def load():
    a = pd.concat([pd.read_csv(p, dtype={"size_filter": str}) for p in sorted((ROOT / "eval").glob("*.csv"))])
    b = pd.concat([pd.read_csv(p, dtype={"size_filter": str}) for p in sorted((ROOT / "eval_v3").glob("*.csv"))])
    if sorted(a.crop_id.unique()) != sorted(b.crop_id.unique()) or b.crop_id.nunique() != 52:
        raise RuntimeError(f"crop sets differ or incomplete: phase1f {a.crop_id.nunique()}, v3 {b.crop_id.nunique()}")
    per = pd.concat([a, b.drop(columns=["size_filter_voxels"])], ignore_index=True)
    i1 = pd.DataFrame([json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "eval").glob("*.json"))])
    i3 = pd.DataFrame([json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "eval_v3").glob("*.json"))])
    info = i1.merge(i3.drop(columns=["specimen", "scan", "test_class", "eval_seconds"]), on="crop_id")
    return per, info


def aggregate(per):
    macro, pooled = [], []
    for scope, d in scopes(per):
        for key, g in d.groupby(KEYS, sort=False):
            base = {"scope": scope, **dict(zip(KEYS, key)), "n_crops": len(g)}
            m, p = dict(base), dict(base) | {"valid_gt": int(g.valid_gt_count.sum()),
                                             "pred_eval": int(g.prediction_count_evaluable.sum())}
            for s in ("025", "050", "075"):
                f1 = g[f"f1_iou{s}"].astype(float)
                m |= {f"f1_{s}_mean": f1.mean(), f"f1_{s}_sd": f1.std(ddof=1), f"f1_{s}_median": f1.median(),
                      f"f1_{s}_min": f1.min(), f"f1_{s}_argmin": g.loc[f1.idxmin(), "crop_id"],
                      f"n_below_090_{s}": int((f1 < FLAG_F1).sum())}
                tp, fp, fn = (int(g[f"{k}_iou{s}"].sum()) for k in ("tp", "fp", "fn"))
                pr = tp / (tp + fp) if tp + fp else 0.0
                rc = tp / (tp + fn) if tp + fn else 0.0
                p |= {f"tp_{s}": tp, f"fp_{s}": fp, f"fn_{s}": fn, f"precision_{s}": pr, f"recall_{s}": rc,
                      f"f1_{s}": 2 * pr * rc / (pr + rc) if pr + rc else 0.0}
            macro.append(m)
            pooled.append(p)
    return pd.DataFrame(macro), pd.DataFrame(pooled)


def v3_stages(crop: str) -> dict:
    """GT-free stage quantities of one v3 run (kept intermediates)."""
    base = ROOT / "v3" / crop
    r = {"crashed": (base / "FAILED.json").exists(), "empty": (base / "EMPTY").exists()}
    sj = base / "00_scale" / "scale.json"
    if sj.exists():
        s = json.loads(sj.read_text(encoding="utf-8"))
        r |= {"r0": s["r0"], "specimen_fraction": s["specimen_fraction"], "fg_fraction": s["fg_fraction_in_specimen"]}
    hits = sorted((base / "pipe" / "03_background").glob("*_bgreject_report.csv"))
    if hits:
        d = pd.read_csv(hits[0])
        r |= {"seed_components": len(d), "bg_keep": int((d.reason == "keep").sum()),
              "bg_external": int((d.reason == "external_background").sum())}
    logf = base / "commands.log"
    if logf.exists():
        txt = logf.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"otsu=(\d+), percentile[\d.]+=([\d.]+), dark_threshold=(\d+)", txt)
        if m:
            r["bg_dark_thr"] = int(m[3])
        for k in ("dark voxels", "external dark", "internal dark"):
            v = re.search(rf"\[bgreject[^\]]*\] {k}\s+=\s+(\d+)", txt)
            r[f"bg_{k.replace(' ', '_')}"] = int(v[1]) if v else None
        if r.get("bg_dark_voxels"):
            r["bg_dark_fraction"] = r["bg_dark_voxels"] / 200 ** 3
            r["bg_external_share"] = r["bg_external_dark"] / r["bg_dark_voxels"]
    s = base / "e11" / "spatial" / "spatialization_summary.json"
    if s.exists():
        r["internal_components"] = json.loads(s.read_text(encoding="utf-8")).get("internal_components")
    diag = base / "e11" / "spatial" / "candidate_diagnostics_no_gt.csv"
    if diag.exists():
        d = pd.read_csv(diag, usecols=["polarity", "point_internal_component", "region_voxels"])
        m = d[d.polarity == "min"]
        r |= {"min_candidates": len(m), "min_in_internal": int((m.point_internal_component > 0).sum())}
    for k, sub in (("tos_markers", "markers_tos_min"), ("edt_otsu_markers", "markers_edt_otsu")):
        p = base / "e11" / sub / "markers.csv"
        r[k] = int(len(pd.read_csv(p))) if p.exists() else None
    return r


def first_stage(r: dict) -> str:
    """Diagnostic rule (documented in REPORT): first pipeline stage whose GT-free quantity is anomalous."""
    if r.get("crashed"):
        return "crash"
    if r.get("empty"):
        return "e11 empty (no MIN candidate in internal_dark)"
    if r.get("r_ratio") is not None and abs(r["r_ratio"] - 1) > SCALE_TOL:
        return "scale (r_granulo vs GT r_in)"
    if r.get("seed_components") and r.get("bg_keep", 0) / r["seed_components"] < 0.05:
        return "stage-03 background (< 5 % seeds kept)"
    if r.get("tos_markers") is not None and r.get("valid_gt") and r["tos_markers"] < 0.5 * r["valid_gt"]:
        return "e11 markers (< 0.5 x valid GT)"
    if r.get("gate") in (True, "True"):
        return "gate fired (EDT-Otsu route)"
    if r.get("precision") is not None and r.get("recall") is not None:
        return "watershed: " + ("under-segmentation / FN" if r["recall"] < r["precision"] else "over-segmentation / FP")
    return "unclassified"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    per, info = load()
    per.to_csv(OUT / "per_crop_all_modes.csv", index=False)
    info.to_csv(OUT / "crop_info.csv", index=False)
    macro, pooled = aggregate(per)
    macro.to_csv(OUT / "crop_macro.csv", index=False)
    pooled.to_csv(OUT / "pooled.csv", index=False)

    prim = per[(per["mode"] == "markers_dominant") & (per.size_filter == "none")]
    wide = prim.pivot_table(index="crop_id", columns=["method", "route"], values="f1_iou050")
    wide.columns = [f"{a}_{b[7:]}_f1_050" for a, b in wide.columns]
    pr = prim[prim.method == "v3"].set_index(["crop_id", "route"])
    wide = wide.join(info.set_index("crop_id"))
    wide.to_csv(OUT / "per_crop_primary.csv")

    # scale: r_granulo vs GT inscribed radius (also e13 and the Phase-1f granulometric diagnostic)
    sc = info[["crop_id", "specimen", "test_class", "gt_inscribed_r_median", "r_est_e13", "r_ps_peak",
               "v3_r_granulo", "v3_openrec_radius", "openrec_radius"]].copy()
    sc["ratio_granulo"] = sc.v3_r_granulo / sc.gt_inscribed_r_median
    sc["ratio_e13"] = sc.r_est_e13 / sc.gt_inscribed_r_median
    sc.to_csv(OUT / "scale.csv", index=False)
    sc_spec = sc.groupby("specimen").agg(n=("crop_id", "size"), gt_r_in=("gt_inscribed_r_median", "median"),
                                         r_granulo=("v3_r_granulo", "median"), r_granulo_min=("v3_r_granulo", "min"),
                                         r_granulo_max=("v3_r_granulo", "max"), r_e13=("r_est_e13", "median"),
                                         ratio_granulo=("ratio_granulo", "median"), ratio_e13=("ratio_e13", "median"))
    ok = sc.dropna(subset=["v3_r_granulo", "gt_inscribed_r_median"])
    corr = float(np.corrcoef(ok.v3_r_granulo, ok.gt_inscribed_r_median)[0, 1])
    corr_e13 = float(np.corrcoef(ok.r_est_e13, ok.gt_inscribed_r_median)[0, 1])

    # runtime (seconds, end to end per crop; contention noted in the report)
    rt = pd.DataFrame({"crop_id": info.crop_id, "specimen": info.specimen, "frozen_s": info.frozen_seconds,
                       "fixflat_s": info.fixflat_seconds_end_to_end, "v3_s": info.v3_seconds,
                       "frozen_peak_mb": info.frozen_peak_mb, "fixflat_peak_mb": info.fixflat_peak_mb,
                       "v3_peak_mb": info.v3_peak_mb, "v3_openrec_impl": info.v3_openrec_impl})
    rt.to_csv(OUT / "runtime.csv", index=False)

    # forensics for flagged crops (and all crops, for completeness)
    frows = []
    for crop in info.crop_id:
        i = info.set_index("crop_id").loc[crop]
        r = {"crop_id": crop, "specimen": i.specimen, "test_class": i.test_class, **v3_stages(crop),
             "r_granulo": i.v3_r_granulo, "gt_r_in": i.gt_inscribed_r_median, "openrec_radius": i.v3_openrec_radius,
             "stage03_thr": i.v3_stage03_thr, "gate": i.v3_gate_triggered, "mv": i.v3_mv, "qv": i.v3_qv,
             "valid_gt": int(pr.loc[(crop, "BA-TOS_sel"), "valid_gt_count"]),
             "precision": float(pr.loc[(crop, "BA-TOS_sel"), "precision_iou050"]),
             "recall": float(pr.loc[(crop, "BA-TOS_sel"), "recall_iou050"]),
             "f1_tos": wide.loc[crop, "v3_ToS_f1_050"], "f1_sel": wide.loc[crop, "v3_sel_f1_050"],
             "fixflat_f1_sel": wide.loc[crop, "fixflat_sel_f1_050"], "frozen_f1_sel": wide.loc[crop, "frozen_sel_f1_050"]}
        r["r_ratio"] = r["r_granulo"] / r["gt_r_in"] if r["r_granulo"] and r["gt_r_in"] else None
        r["flag"] = bool(min(r["f1_tos"], r["f1_sel"]) < FLAG_F1)
        r["first_stage"] = first_stage(r) if r["flag"] else ""
        frows.append(r)
    fr = pd.DataFrame(frows)
    fr.to_csv(OUT / "forensics_stages.csv", index=False)
    flags = fr[fr.flag]
    flags.to_csv(OUT / "flags.csv", index=False)

    # ---------------- markdown tables
    L = []
    order = {"frozen": 0, "fixflat": 1, "v3": 2}
    for flt in ("none", "500", "own"):
        for kind, tab in (("Crop-macro F1", macro), ("Pooled objects", pooled)):
            t = tab[(tab["mode"] == "markers_dominant") & (tab.size_filter == flt)].copy()
            if t.empty:
                continue
            t["o"] = t.method.map(order)
            t = t.sort_values(["scope", "route", "o"], kind="stable")
            note = " (v3 only: 500 s^3)" if flt == "own" else ""
            L.append(f"\n### {kind} (markers/dominant = seed-ignore), size filter {flt}{note}\n")
            if kind.startswith("Crop"):
                L.append("| scope | route | method | n | F1@0.50 mean ± SD | median | min (crop) | # < 0.90 | F1@0.75 | F1@0.25 |")
                L.append("|---|---|---|---|---|---|---|---|---|---|")
                for r in t.itertuples():
                    L.append(f"| {r.scope} | {r.route} | {r.method} | {r.n_crops} | {r.f1_050_mean:.3f} ± {r.f1_050_sd:.3f} "
                             f"| {r.f1_050_median:.3f} | {r.f1_050_min:.3f} ({r.f1_050_argmin}) | {r.n_below_090_050} "
                             f"| {r.f1_075_mean:.3f} | {r.f1_025_mean:.3f} |")
            else:
                L.append("| scope | route | method | GT | pred | TP | FP | FN | P@0.50 | R@0.50 | F1@0.50 | F1@0.75 |")
                L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
                for r in t.itertuples():
                    L.append(f"| {r.scope} | {r.route} | {r.method} | {r.valid_gt} | {r.pred_eval} | {r.tp_050} | {r.fp_050} "
                             f"| {r.fn_050} | {r.precision_050:.3f} | {r.recall_050:.3f} | {r.f1_050:.3f} | {r.f1_075:.3f} |")
    t = macro[(macro["mode"] == "no_markers_fraction050") & macro.scope.isin(["all", "EFR_internal", "external"])].copy()
    t["o"] = t.method.map(order)
    L.append("\n### No-marker mode (own-voxel support, overlap-ignore 0.50), crop-macro F1@0.50\n")
    L.append("| scope | route | filter | method | mean ± SD | min |")
    L.append("|---|---|---|---|---|---|")
    for r in t.sort_values(["scope", "route", "size_filter", "o"]).itertuples():
        L.append(f"| {r.scope} | {r.route} | {r.size_filter} | {r.method} | {r.f1_050_mean:.3f} ± {r.f1_050_sd:.3f} "
                 f"| {r.f1_050_min:.3f} |")

    L.append("\n### Gate (E23c; v3 limits 3·max(1,s³) / 4·max(1,s³))\n")
    L.append("| specimen | n | frozen fires | fixflat fires | v3 fires | frozen crash | fixflat crash | v3 crash / empty |")
    L.append("|---|---|---|---|---|---|---|---|")
    for spec, g in list(info.groupby("specimen")) + [("all", info)]:
        fire = lambda c: int((g[c].astype(str) == "True").sum())  # noqa: E731
        L.append(f"| {spec} | {len(g)} | {fire('frozen_gate_triggered')} | {fire('fixflat_gate_triggered')} "
                 f"| {fire('v3_gate_triggered')} | {int(g.frozen_failed.notna().sum())} | {int(g.fixflat_failed.notna().sum())} "
                 f"| {int(g.v3_failed.notna().sum())} / {int(g.v3_empty.notna().sum())} |")

    L.append(f"\n### Scale: r_granulo vs GT median inscribed radius (GT diagnostic only); Pearson r = {corr:.3f} "
             f"(e13: {corr_e13:.3f}), n = {len(ok)}\n")
    L.append("| specimen | n | GT r_in (median) | r_granulo median [min–max] | r_granulo / GT r_in | r_e13 median | r_e13 / GT r_in |")
    L.append("|---|---|---|---|---|---|---|")
    for spec, r in sc_spec.iterrows():
        L.append(f"| {spec} | {r.n} | {r.gt_r_in:.2f} | {r.r_granulo:g} [{r.r_granulo_min:g}–{r.r_granulo_max:g}] "
                 f"| {r.ratio_granulo:.2f} | {r.r_e13:.2f} | {r.ratio_e13:.2f} |")

    L.append("\n### Runtime per crop (s, end to end; peak MB of the largest stage)\n")
    L.append("| specimen | n | frozen s (median) | fixflat s | v3 s | v3 s max | v3 peak MB max | v3 lowmem openrec |")
    L.append("|---|---|---|---|---|---|---|---|")
    for spec, g in list(rt.groupby("specimen")) + [("all", rt)]:
        L.append(f"| {spec} | {len(g)} | {g.frozen_s.median():.0f} | {g.fixflat_s.median():.0f} | {g.v3_s.median():.0f} "
                 f"| {g.v3_s.max():.0f} | {g.v3_peak_mb.max():.0f} | {int((g.v3_openrec_impl == 'openrec_lowmem').sum())} |")
    L.append(f"\nTotals (h): frozen {rt.frozen_s.sum() / 3600:.1f}, fixflat {rt.fixflat_s.sum() / 3600:.1f}, "
             f"v3 {rt.v3_s.sum() / 3600:.1f}\n")

    L.append(f"\n### Flagged crops: v3 F1@0.50 < {FLAG_F1} (ToS or sel), first-stage forensics\n")
    L.append("| crop | class | v3 ToS / sel | fixflat sel | frozen sel | r_granulo / GT r_in | spec. frac | openrec r "
             "| thr03 | seeds kept / total | ext. dark share | MIN in internal | ToS markers / valid GT | gate (mv/qv) "
             "| P / R | first stage |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    fmt = lambda v, f="{:g}": "—" if v is None or (isinstance(v, float) and np.isnan(v)) else f.format(v)  # noqa: E731
    for r in flags.to_dict("records"):
        L.append(f"| {r['crop_id']} | {r['test_class']} | {r['f1_tos']:.3f} / {r['f1_sel']:.3f} | {r['fixflat_f1_sel']:.3f} "
                 f"| {r['frozen_f1_sel']:.3f} | {fmt(r['r_granulo'])} / {fmt(r['gt_r_in'], '{:.1f}')} "
                 f"| {fmt(r.get('specimen_fraction'), '{:.2f}')} | {fmt(r['openrec_radius'])} | {fmt(r['stage03_thr'])} "
                 f"| {fmt(r.get('bg_keep'))} / {fmt(r.get('seed_components'))} | {fmt(r.get('bg_external_share'), '{:.2f}')} "
                 f"| {fmt(r.get('min_in_internal'))} | {fmt(r.get('tos_markers'))} / {r['valid_gt']} "
                 f"| {r['gate']} ({fmt(r['mv'])}/{fmt(r['qv'])}) | {r['precision']:.2f} / {r['recall']:.2f} "
                 f"| {r['first_stage']} |")

    L.append("\n### Per crop (markers/dominant, no filter)\n")
    L.append("| crop | class | frozen ToS / sel | fixflat ToS / sel | v3 ToS / sel | r e13 / r_granulo / GT r_in | v3 thr03 "
             "| v3 gate (mv/qv) | v3 markers | v3 s |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for cid, r in wide.iterrows():
        st = "CRASH" if isinstance(r.v3_failed, str) else ("EMPTY" if isinstance(r.v3_empty, str) else
                                                            f"{'FIRES' if str(r.v3_gate_triggered) == 'True' else 'off'} ({r.v3_mv:g}/{r.v3_qv:g})")
        L.append(f"| {cid} | {r.test_class} | {r.frozen_ToS_f1_050:.3f} / {r.frozen_sel_f1_050:.3f} | {r.fixflat_ToS_f1_050:.3f}"
                 f" / {r.fixflat_sel_f1_050:.3f} | {r.v3_ToS_f1_050:.3f} / {r.v3_sel_f1_050:.3f} | {r.r_est_e13:.1f} / "
                 f"{fmt(r.v3_r_granulo)} / {r.gt_inscribed_r_median:.1f} | {fmt(r.v3_stage03_thr)} | {st} "
                 f"| {fmt(r.v3_tos_markers)} | {fmt(r.v3_seconds, '{:.0f}')} |")
    (OUT / "tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    git = git_state()["commit"]
    rows = []
    ii = info.set_index("crop_id")
    for r in per[per.method == "v3"].itertuples():
        i = ii.loc[r.crop_id]
        rows.append({"experiment_id": f"test_v3_v3/{r.crop_id}", "method": r.route, "model": "batos-paper-v3-scalefree",
                     "volume": r.scan, "crop": r.crop_id, "device": "cpu (AMD Ryzen 7 5700X)",
                     "parameters": f"size_filter={r.size_filter};eval=evaluate_v2:{r.mode};split=test;"
                                   f"openrec={i.v3_openrec_impl};r_granulo={i.v3_r_granulo}",
                     "runtime_seconds": i.v3_seconds, "peak_memory_mb": i.v3_peak_mb, "tp": r.tp_iou050,
                     "fp": r.fp_iou050, "fn": r.fn_iou050, "precision": r.precision_iou050, "recall": r.recall_iou050,
                     "f1_025": r.f1_iou025, "f1_050": r.f1_iou050, "f1_075": r.f1_iou075, "ap50": r.ap_iou050,
                     "ap75": r.ap_iou075, "mean_matched_iou": r.mean_matched_iou050,
                     "n_pred": r.prediction_count_evaluable, "n_gt": r.valid_gt_count, "git_commit": git,
                     "run_json": str(ROOT / "eval_v3" / f"{r.crop_id}.json")})
    upsert(REPO / "results" / "all_runs.csv", rows, ["experiment_id", "method", "parameters"])
    print((OUT / "tables.md").read_text(encoding="utf-8"))
    print(f"flags (v3 F1@0.50 < {FLAG_F1}): {list(flags.crop_id)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
