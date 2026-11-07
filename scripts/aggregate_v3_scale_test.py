#!/usr/bin/env python3
"""Phase 5e step 3: aggregate the synthetic scale test (scripts/run_v3_scale_test.py) into results/v3_scale_test/.

Crop-macro (per factor x method x route x size filter x mode) and pooled objects are kept in separate tables.
Outputs: per_crop_all_modes.csv, crop_info.csv, crop_macro.csv, pooled.csv, per_crop_primary.csv, tables.md.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
EVAL = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test\eval")
OUT = REPO / "results" / "v3_scale_test"
KEYS = ["factor", "method", "route", "size_filter", "mode"]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    per = pd.concat([pd.read_csv(p) for p in sorted(EVAL.glob("*.csv"))], ignore_index=True)
    info = pd.DataFrame([json.loads(p.read_text(encoding="utf-8")) for p in sorted(EVAL.glob("*.json"))])
    per.to_csv(OUT / "per_crop_all_modes.csv", index=False)
    info.to_csv(OUT / "crop_info.csv", index=False)
    macro, pooled = [], []
    for key, g in per.groupby(KEYS, sort=True):
        base = {**dict(zip(KEYS, key)), "n_crops": len(g)}
        m, p = dict(base), dict(base) | {"valid_gt": int(g.valid_gt_count.sum()),
                                         "pred_eval": int(g.prediction_count_evaluable.sum())}
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
    wide = prim.pivot_table(index=["dev_crop_id", "factor"], columns=["method", "route"], values="f1_iou050")
    wide.columns = [f"{a}_{b[7:]}_f1_050" for a, b in wide.columns]
    wide = wide.reset_index().merge(info, on=["dev_crop_id", "factor"], how="left")
    wide.to_csv(OUT / "per_crop_primary.csv", index=False)

    lines = []
    for flt in ("none", "own"):
        t = macro[(macro["mode"] == "markers_dominant") & (macro.size_filter == flt)]
        lines.append(f"\n### Crop-macro F1 (markers/dominant = seed-ignore), size filter {flt}"
                     f"{' (fixflat 500, v3 500 s^3)' if flt == 'own' else ''}\n")
        lines.append("| factor | method | route | n | F1@0.50 mean ± SD | median | min (crop) | F1@0.75 mean | F1@0.25 mean |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for r in t.itertuples():
            lines.append(f"| {r.factor:g} | {r.method} | {r.route} | {r.n_crops} | {r.f1_050_mean:.3f} ± {r.f1_050_sd:.3f}"
                         f" | {r.f1_050_median:.3f} | {r.f1_050_min:.3f} ({r.f1_050_argmin}) | {r.f1_075_mean:.3f} "
                         f"| {r.f1_025_mean:.3f} |")
        t = pooled[(pooled["mode"] == "markers_dominant") & (pooled.size_filter == flt)]
        lines.append(f"\n### Pooled objects, size filter {flt}\n")
        lines.append("| factor | method | route | GT | pred | TP | FP | FN | P@0.50 | R@0.50 | F1@0.50 | F1@0.75 |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in t.itertuples():
            lines.append(f"| {r.factor:g} | {r.method} | {r.route} | {r.valid_gt} | {r.pred_eval} | {r.tp_050} | {r.fp_050} "
                         f"| {r.fn_050} | {r.precision_050:.3f} | {r.recall_050:.3f} | {r.f1_050:.3f} | {r.f1_075:.3f} |")
    t = macro[macro["mode"] == "no_markers_fraction050"]
    lines.append("\n### No-marker mode (own-voxel support, overlap-ignore 0.50), crop-macro F1@0.50\n")
    lines.append("| factor | method | route | filter | mean ± SD | min |")
    lines.append("|---|---|---|---|---|---|")
    for r in t.itertuples():
        lines.append(f"| {r.factor:g} | {r.method} | {r.route} | {r.size_filter} | {r.f1_050_mean:.3f} ± {r.f1_050_sd:.3f} "
                     f"| {r.f1_050_min:.3f} |")
    lines.append("\n### Per crop (markers/dominant, no filter)\n")
    lines.append("| crop | f | fixflat ToS / sel | v3 ToS / sel | r e13 / r_granulo | GT r_in | openrec fixflat / v3 "
                 "| thr03 v3 | gate fixflat (mv/qv) | gate v3 (mv/qv) | markers fixflat / v3 | v3 s / peak MB |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")

    def gate(r, m):
        if isinstance(r.get(f"{m}_failed"), str):
            return "CRASH"
        if isinstance(r.get(f"{m}_empty"), str):
            return "EMPTY"
        t = r.get(f"{m}_gate")
        return f"{'FIRES' if t in (True, 'True') else 'off'} ({r.get(f'{m}_mv')}/{r.get(f'{m}_qv')})"

    for r in wide.sort_values(["factor", "dev_crop_id"]).to_dict("records"):
        lines.append(f"| {r['dev_crop_id']} | {r['factor']:g} | {r['fixflat_ToS_f1_050']:.3f} / {r['fixflat_sel_f1_050']:.3f}"
                     f" | {r['v3_ToS_f1_050']:.3f} / {r['v3_sel_f1_050']:.3f} | {r['fixflat_r_scale']} / {r['v3_r_scale']}"
                     f" | {r['gt_inscribed_r_median']:.2f} | {r['fixflat_openrec_radius']} / {r['v3_openrec_radius']} "
                     f"| {r['v3_stage03_thr']} | {gate(r, 'fixflat')} | {gate(r, 'v3')} | {r['fixflat_tos_markers']} / "
                     f"{r['v3_tos_markers']} | {r['v3_seconds']} / {r['v3_peak_mb']} |")
    (OUT / "tables.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print((OUT / "tables.md").read_text(encoding="utf-8")[:6000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
