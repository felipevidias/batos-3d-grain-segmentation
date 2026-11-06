#!/usr/bin/env python3
"""Phase 1f: aggregate the 52 test crops (scripts/run_test_protocol.py outputs) into results/test_protocol/.

Inputs : C:\\BA_TOS_DATA\\runs\\test_protocol\\eval\\<crop>.csv (evaluate_v2 rows) and <crop>.json (scale, gate, runtime).
Outputs: per_crop_all_modes.csv, crop_info.csv, crop_macro.csv, pooled.csv, flags.csv, tables.md;
         rows in results/all_runs.csv (experiment_id test_protocol_<method>/<crop>; key experiment_id+method+parameters).
Scopes : each specimen, EFR_internal (= EFR), external (5 specimens), all (52). Crop-macro and pooled kept separate.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.extinction_run import upsert  # noqa: E402
from batos_v2.run_crop import git_state  # noqa: E402

EVAL = Path(r"C:\BA_TOS_DATA\runs\test_protocol\eval")
OUT = REPO / "results" / "test_protocol"
KEYS = ["method", "route", "size_filter", "mode"]
FLAG_F1 = 0.90


def scopes(df: pd.DataFrame):
    yield "all", df
    yield "EFR_internal", df[df.specimen == "EFR"]
    yield "external", df[df.specimen != "EFR"]
    for s in sorted(df.specimen.unique()):
        yield s, df[df.specimen == s]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    per = pd.concat([pd.read_csv(p, dtype={"size_filter": str}) for p in sorted(EVAL.glob("*.csv"))],
                    ignore_index=True)
    info = pd.DataFrame([json.loads(p.read_text(encoding="utf-8")) for p in sorted(EVAL.glob("*.json"))])
    per.to_csv(OUT / "per_crop_all_modes.csv", index=False)
    info.to_csv(OUT / "crop_info.csv", index=False)
    macro, pooled = [], []
    for scope, d in scopes(per):
        for key, g in d.groupby(KEYS, sort=False):
            base = {"scope": scope, **dict(zip(KEYS, key)), "n_crops": len(g)}
            m, p = dict(base), dict(base)
            p |= {"valid_gt": int(g.valid_gt_count.sum()), "pred_eval": int(g.prediction_count_evaluable.sum())}
            for s in ("025", "050", "075"):
                f1 = g[f"f1_iou{s}"].astype(float)
                m |= {f"f1_{s}_mean": f1.mean(), f"f1_{s}_sd": f1.std(ddof=1), f"f1_{s}_median": f1.median(),
                      f"f1_{s}_min": f1.min(), f"f1_{s}_argmin": g.loc[f1.idxmin(), "crop_id"]}
                tp, fp, fn = (int(g[f"{k}_iou{s}"].sum()) for k in ("tp", "fp", "fn"))
                pr = tp / (tp + fp) if tp + fp else 0.0
                rc = tp / (tp + fn) if tp + fn else 0.0
                p |= {f"tp_{s}": tp, f"fp_{s}": fp, f"fn_{s}": fn, f"precision_{s}": pr, f"recall_{s}": rc,
                      f"f1_{s}": 2 * pr * rc / (pr + rc) if pr + rc else 0.0}
            macro.append(m)
            pooled.append(p)
    macro, pooled = pd.DataFrame(macro), pd.DataFrame(pooled)
    macro.to_csv(OUT / "crop_macro.csv", index=False)
    pooled.to_csv(OUT / "pooled.csv", index=False)

    prim = per[(per["mode"] == "markers_dominant") & (per.size_filter == "none")]
    wide = prim.pivot_table(index="crop_id", columns=["method", "route"], values="f1_iou050")
    wide.columns = [f"{a}_{b[7:]}_f1_050" for a, b in wide.columns]
    wide = wide.join(info.set_index("crop_id"))
    flags = wide[(wide["fixflat_ToS_f1_050"] < FLAG_F1) | (wide["fixflat_sel_f1_050"] < FLAG_F1)]
    flags.to_csv(OUT / "flags.csv")
    wide.to_csv(OUT / "per_crop_primary.csv")

    # markdown tables (primary mode = markers_dominant)
    lines = []
    for flt in ("none", "500"):
        t = macro[(macro["mode"] == "markers_dominant") & (macro.size_filter == flt)]
        lines.append(f"\n### Crop-macro F1 (evaluate_v2 markers/dominant = seed-ignore), size filter {flt}\n")
        lines.append("| scope | method | route | n | F1@0.50 mean ± SD | median | min (crop) | F1@0.75 mean | F1@0.25 mean |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for r in t.itertuples():
            lines.append(f"| {r.scope} | {r.method} | {r.route} | {r.n_crops} | {r.f1_050_mean:.3f} ± {r.f1_050_sd:.3f} "
                         f"| {r.f1_050_median:.3f} | {r.f1_050_min:.3f} ({r.f1_050_argmin}) | {r.f1_075_mean:.3f} "
                         f"| {r.f1_025_mean:.3f} |")
        t = pooled[(pooled["mode"] == "markers_dominant") & (pooled.size_filter == flt)]
        lines.append(f"\n### Pooled objects, size filter {flt}\n")
        lines.append("| scope | method | route | GT | pred | TP | FP | FN | P@0.50 | R@0.50 | F1@0.50 | F1@0.75 |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in t.itertuples():
            lines.append(f"| {r.scope} | {r.method} | {r.route} | {r.valid_gt} | {r.pred_eval} | {r.tp_050} | {r.fp_050} "
                         f"| {r.fn_050} | {r.precision_050:.3f} | {r.recall_050:.3f} | {r.f1_050:.3f} | {r.f1_075:.3f} |")
    t = macro[(macro["mode"] == "no_markers_fraction050") & (macro.scope.isin(["all", "EFR_internal", "external"]))]
    lines.append("\n### No-marker mode (own-voxel support, overlap-ignore 0.50), crop-macro F1@0.50\n")
    lines.append("| scope | method | route | filter | mean ± SD | min |")
    lines.append("|---|---|---|---|---|---|")
    for r in t.itertuples():
        lines.append(f"| {r.scope} | {r.method} | {r.route} | {r.size_filter} | {r.f1_050_mean:.3f} ± {r.f1_050_sd:.3f} "
                     f"| {r.f1_050_min:.3f} |")
    lines.append("\n### Per crop (markers/dominant, no filter), gate, scale, runtime\n")
    lines.append("| crop | class | f | frozen ToS / sel F1@0.50 | fixflat ToS / sel | gate frozen (mv/qv) | gate fixflat (mv/qv) "
                 "| r_est e13 | r_ps_peak | GT r_in | openrec r (impl) | frozen s / peak MB | fixflat s / peak MB |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for cid, r in wide.iterrows():
        gf = f"{'FIRES' if r.frozen_gate_triggered else 'off'} ({r.frozen_mv:g}/{r.frozen_qv:g})"
        gx = f"{'FIRES' if r.fixflat_gate_triggered else 'off'} ({r.fixflat_mv:g}/{r.fixflat_qv:g})"
        impl = "lowmem" if r.openrec_implementation == "openrec_lowmem" else "frozen"
        lines.append(f"| {cid} | {r.test_class} | {r.specimen_fraction_sub4:.2f} | {r.frozen_ToS_f1_050:.3f} / "
                     f"{r.frozen_sel_f1_050:.3f} | {r.fixflat_ToS_f1_050:.3f} / {r.fixflat_sel_f1_050:.3f} | {gf} | {gx} "
                     f"| {r.r_est_e13:.2f} | {r.r_ps_peak} | {r.gt_inscribed_r_median:.2f} | {r.openrec_radius} ({impl}) "
                     f"| {r.frozen_seconds:.0f} / {r.frozen_peak_mb:.0f} | {r.fixflat_seconds_end_to_end:.0f} / "
                     f"{r.fixflat_peak_mb:.0f} |")
    (OUT / "tables.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    git = git_state()["commit"]
    rows = []
    for r in per.itertuples():
        i = info.set_index("crop_id").loc[r.crop_id]
        secs = i.frozen_seconds if r.method == "frozen" else i.fixflat_seconds_end_to_end
        peak = i.frozen_peak_mb if r.method == "frozen" else i.fixflat_peak_mb
        rows.append({"experiment_id": f"test_protocol_{r.method}/{r.crop_id}", "method": r.route,
                     "model": "batos-paper-v2-frozen" if r.method == "frozen" else "batos-paper-v2-fixflat",
                     "volume": r.scan, "crop": r.crop_id, "device": "cpu (AMD Ryzen 7 5700X)",
                     "parameters": f"size_filter={r.size_filter};eval=evaluate_v2:{r.mode};split=test;"
                                   f"openrec={i.openrec_implementation}",
                     "runtime_seconds": secs, "peak_memory_mb": peak, "tp": r.tp_iou050, "fp": r.fp_iou050,
                     "fn": r.fn_iou050, "precision": r.precision_iou050, "recall": r.recall_iou050,
                     "f1_025": r.f1_iou025, "f1_050": r.f1_iou050, "f1_075": r.f1_iou075, "ap50": r.ap_iou050,
                     "ap75": r.ap_iou075, "mean_matched_iou": r.mean_matched_iou050,
                     "n_pred": r.prediction_count_evaluable, "n_gt": r.valid_gt_count, "git_commit": git,
                     "run_json": str(EVAL / f"{r.crop_id}.json")})
    upsert(REPO / "results" / "all_runs.csv", rows, ["experiment_id", "method", "parameters"])
    print((OUT / "tables.md").read_text(encoding="utf-8"))
    print(f"flags (fixflat F1@0.50 < {FLAG_F1}): {list(flags.index)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
