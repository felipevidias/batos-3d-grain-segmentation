#!/usr/bin/env python3
"""Phase 5h: aggregate the v3.2 uniqueness hypotheses on dev + synthetic crops against v3 and v3.1.

v3, fixflat: Phase 5e step-3 evaluations (C:\\BA_TOS_DATA\\runs\\v3_scale_test\\eval). v3.1: Phase 5f
(C:\\BA_TOS_DATA\\runs\\v31_dev\\eval for f = 1.4 / 2.0; = v3 for f <= 1). v3.2-<H>: C:\\BA_TOS_DATA\\runs\\v32_dev\\eval\\<H>.
Seed diagnostics (markers per grain, seeded, split) for v3.2 come from its evaluations; for v3 / v3.1 they are recomputed
here from the kept normalized markers (GT = diagnosis only). Acceptance per docs/V32_UNIQUENESS.md.
Outputs: results/v32_dev/{crop_macro,pooled,per_crop,seed_density,acceptance}.csv, tables.md.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import tifffile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from run_v32_dev import seed_stats  # noqa: E402
from run_v3_scale_test import crop_sets  # noqa: E402

V3 = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test")
V31 = Path(r"C:\BA_TOS_DATA\runs\v31_dev")
V32 = Path(r"C:\BA_TOS_DATA\runs\v32_dev")
OUT = REPO / "results" / "v32_dev"
V3_REF = {0.7: 0.975, 1.0: 0.989, 1.4: 0.995, 2.0: 0.938}  # docs/V32_UNIQUENESS.md (v3 crop-macro, rounded as reported)
TOL = 0.005
SEL = "e23c/watershed/normalized_markers/batos_selective.tif"


def load(hyps: list[str]) -> pd.DataFrame:
    rd = lambda ps: pd.concat([pd.read_csv(p, dtype={"size_filter": str}) for p in ps])  # noqa: E731
    base = rd(sorted((V3 / "eval").glob("*.csv")))
    v31 = rd(sorted((V31 / "eval").glob("*.csv")))
    low = base[(base.method == "v3") & (base.factor <= 1.0)].assign(method="v3.1")
    parts = [base, low, v31]
    for h in hyps:
        ps = sorted((V32 / "eval" / h).glob("*.csv"))
        if ps:
            parts.append(rd(ps))
    return pd.concat(parts, ignore_index=True)


def ref_seeds(method: str, f: float, c: dict) -> dict:
    root = V3 / "runs" / "v3" if (method == "v3" or f <= 1.0) else V31 / "runs"
    p = root / f"x{f:g}" / c["crop_id"] / SEL
    if not p.exists():
        return {}
    return seed_stats(tifffile.imread(p), tifffile.imread(c["gt_path"]), c["min_gt_voxels"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--hyps", nargs="*", default=["H1", "H2", "H3"])
    a = ap.parse_args(argv)
    per = load(a.hyps)
    methods = [m for m in ["v3", "v3.1"] + [f"v3.2-{h}" for h in a.hyps] if m in set(per.method)]
    OUT.mkdir(parents=True, exist_ok=True)
    keys = ["factor", "method", "route", "size_filter", "mode"]
    macro, pooled = [], []
    for k, g in per[per.method.isin(methods)].groupby(keys):
        b = dict(zip(keys, k)) | {"n_crops": len(g)}
        m = dict(b)
        for s in ("025", "050", "075"):
            f1 = g[f"f1_iou{s}"]
            m |= {f"f1_{s}_mean": f1.mean(), f"f1_{s}_sd": f1.std(ddof=1), f"f1_{s}_min": f1.min(),
                  f"f1_{s}_argmin": g.loc[f1.idxmin(), "crop_id"], f"n_below_090_{s}": int((f1 < 0.90).sum())}
        macro.append(m)
        row = dict(b)
        for s in ("050", "075"):
            tp, fp, fn = (int(g[f"{x}_iou{s}"].sum()) for x in ("tp", "fp", "fn"))
            row |= {f"tp_{s}": tp, f"fp_{s}": fp, f"fn_{s}": fn, f"f1_{s}": 2 * tp / (2 * tp + fp + fn) if tp else 0.0}
        pooled.append(row | {"valid_gt": int(g.valid_gt_count.sum())})
    macro, pooled = pd.DataFrame(macro), pd.DataFrame(pooled)
    macro.to_csv(OUT / "crop_macro.csv", index=False)
    pooled.to_csv(OUT / "pooled.csv", index=False)

    prim = per[(per["mode"] == "markers_dominant") & (per.size_filter == "none") & (per.route == "BA-TOS_sel")
               & per.method.isin(methods)]
    sets = crop_sets()
    cmap = {c["crop_id"]: c for f in sets for c in sets[f]}
    rows = []
    for r in prim.itertuples():
        if r.method.startswith("v3.2"):
            sd = {k[5:]: getattr(r, k) for k in r._fields if k.startswith("seed_")}
        else:
            sd = ref_seeds(r.method, r.factor, cmap[r.crop_id])
        rows.append({"crop_id": r.crop_id, "dev_crop_id": r.dev_crop_id, "factor": r.factor, "method": r.method,
                     "f1_050": r.f1_iou050, "f1_075": r.f1_iou075, "precision": r.precision_iou050,
                     "recall": r.recall_iou050, "pred_total": r.n_pred_total, **sd})
    sd = pd.DataFrame(rows)
    sd.to_csv(OUT / "seed_density.csv", index=False)
    wide = prim.pivot_table(index=["dev_crop_id", "factor"], columns="method", values="f1_iou050").reset_index()
    wide.to_csv(OUT / "per_crop.csv", index=False)

    acc = []
    for m in [x for x in methods if x.startswith("v3.2")]:
        for f in sorted(V3_REF):
            q = macro[(macro.factor == f) & (macro.method == m) & (macro.route == "BA-TOS_sel")
                      & (macro.size_filter == "none") & (macro["mode"] == "markers_dominant")]
            if q.empty:
                acc.append({"method": m, "factor": f, "n_crops": 0, "f1_050": None, "threshold": None, "pass": None})
                continue
            v = float(q.f1_050_mean.iloc[0])
            thr = 0.985 if f == 1.0 else V3_REF[f] - TOL
            acc.append({"method": m, "factor": f, "n_crops": int(q.n_crops.iloc[0]), "f1_050": v, "threshold": thr,
                        "pass": bool(v >= thr) if int(q.n_crops.iloc[0]) == 12 else None})
    acc = pd.DataFrame(acc)
    acc.to_csv(OUT / "acceptance.csv", index=False)

    L = ["# v3.2 uniqueness hypotheses on dev + synthetic crops (generated by scripts/aggregate_v32_dev.py)", "",
         "Primary: crop-macro F1 BA-TOS_sel, seed-ignore (markers_dominant). v3.1 = v3 for f <= 1.", ""]
    L += ["## Acceptance (docs/V32_UNIQUENESS.md)", "", "| method | f | n | F1@0.50 | threshold | pass |",
          "|---|---|---|---|---|---|"]
    for r in acc.itertuples():
        L.append(f"| {r.method} | {r.factor:g} | {r.n_crops} | {'' if r.f1_050 is None else f'{r.f1_050:.3f}'} | "
                 f"{'' if r.threshold is None else f'{r.threshold:.3f}'} | {r._6} |")
    for flt in ("none", "own"):
        for s in ("050", "075"):
            L += ["", f"## Crop-macro F1@0.{s[1:]}, size filter {flt} (mean ± SD, min, # < 0.90; pooled)", "",
                  "| f | " + " | ".join(methods) + " |", "|---|" + "---|" * len(methods)]
            for f in sorted(per.factor.unique()):
                cells = []
                for m in methods:
                    q = macro[(macro.factor == f) & (macro.method == m) & (macro.route == "BA-TOS_sel")
                              & (macro.size_filter == flt) & (macro["mode"] == "markers_dominant")]
                    p = pooled[(pooled.factor == f) & (pooled.method == m) & (pooled.route == "BA-TOS_sel")
                               & (pooled.size_filter == flt) & (pooled["mode"] == "markers_dominant")]
                    cells.append("—" if q.empty else
                                 f"{q[f'f1_{s}_mean'].iloc[0]:.3f} ± {q[f'f1_{s}_sd'].iloc[0]:.3f}, min "
                                 f"{q[f'f1_{s}_min'].iloc[0]:.3f}, {int(q[f'n_below_090_{s}'].iloc[0])}; "
                                 f"pooled {p[f'f1_{s}'].iloc[0]:.3f} (n {int(q.n_crops.iloc[0])})")
                L.append(f"| {f:g} | " + " | ".join(cells) + " |")
    L += ["", "## Seed diagnostics, BA-TOS_sel markers (medians over crops; GT = diagnosis only)", "",
          "| f | method | markers / grain (all GT) | markers / valid grain | valid seeded | valid split | "
          "markers in background | pred / valid GT |", "|---|---|---|---|---|---|---|---|"]
    for (f, m), g in sd.groupby(["factor", "method"]):
        if "frac_seeded" not in g or g.frac_seeded.isna().all():
            continue
        L.append(f"| {f:g} | {m} | {g.markers_per_grain_all.median():.2f} | {g.markers_per_valid_grain.median():.2f} | "
                 f"{g.frac_seeded.median():.3f} | {g.frac_split.median():.3f} | {g.markers_background.median():.0f} | "
                 f"{(g.pred_total / g.gt_valid).median():.2f} |")
    L += ["", "## Per crop F1@0.50 (sel, no filter)", "",
          "| dev crop | f | " + " | ".join(methods) + " |", "|---|---|" + "---|" * len(methods)]
    for r in wide.sort_values(["factor", "dev_crop_id"]).itertuples(index=False):
        vals = dict(zip(wide.columns, r))
        L.append(f"| {vals['dev_crop_id']} | {vals['factor']:g} | "
                 + " | ".join("—" if pd.isna(vals.get(m)) else f"{vals[m]:.3f}" for m in methods) + " |")
    (OUT / "tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L[:12]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
