#!/usr/bin/env python3
"""Phase 1d verification: evaluate_v2 (geometric seed-ignore) vs the frozen e09b on the 12 frozen runs.

Reference = the Phase 1 closure e09b outputs (frozen e09b run unmodified + IoU 0.25 via
scripts/e09b_with_iou025.py), C:\\BA_TOS_DATA\\runs\\batos_paper_v2_frozen\\<crop>\\phase1_closure\\<label>_<filter>.
For every crop x route (BA-TOS_ToS, BA-TOS_sel, Otsu+WS) x size filter (none, 500):
  A. markers supplied, rule=dominant          -> must equal e09b seed_ignore (all counts and metrics, IoU .25/.50/.75)
  B. same, prediction labels randomly permuted -> must equal A (label invariance; e09b cannot do this)
  C. no markers, rule=fraction, threshold 0.50 -> compared with e09b overlap_ignore (expected identical)
  D. markers, rule=fraction, threshold 0.50   -> reported only
Writes results/evaluate_v2/{verification.csv, per_crop_all_modes.csv,
crop_macro_by_mode.csv, pooled_object_by_mode.csv,
morphology_summary_pooled.csv, morphology_matched_seed.csv.gz (untracked), REPORT.md}. Fails (exit 1) if A or B mismatch.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.evaluate_v2 import evaluate  # noqa: E402

RUNS = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen")
CROPS_DIR = Path(r"C:\BA_TOS_DATA\crops")
OUT = REPO / "results" / "evaluate_v2"
CROPS = [f"EFRGP01_0{s}_crop0{c}" for s in (0, 1) for c in (2, 3, 4, 5, 6, 7)]
ROUTES = {  # route -> (watershed dir, configuration, marker label)
    "BA-TOS_ToS": ("e11/watershed", "tos_e11_min_distance_c1", "tos_e11_min"),
    "BA-TOS_sel": ("e23c/watershed", "batos_selective_distance_c1", "batos_selective"),
    "Otsu+WS(edt_otsu)": ("e11/watershed", "edt_otsu_distance_c1", "edt_otsu"),
}
METRIC_COLS = ["prediction_count_evaluable", "valid_gt_count"] + [
    f"{m}{s}" for s in ("025", "050", "075")
    for m in ("tp_iou", "fp_iou", "fn_iou", "precision_iou", "recall_iou", "f1_iou", "mean_matched_iou")]
INT_COLS = {c for c in METRIC_COLS if c.split("_iou")[0] in ("tp", "fp", "fn") or c.endswith("count")
            or c == "prediction_count_evaluable"}


def same(a: dict, b: pd.Series, skip: tuple[str, ...] = ()) -> tuple[bool, float]:
    worst = 0.0
    ok = True
    for c in METRIC_COLS:
        if c in skip:
            continue
        x, y = float(a[c]), float(b[c])
        if c in INT_COLS:
            ok &= x == y
        else:
            worst = max(worst, abs(x - y))
    return bool(ok and worst <= 1e-12), worst


def md(frame: pd.DataFrame) -> str:
    cells = [[f"{v:.4f}" if isinstance(v, float) else str(v) for v in r] for r in frame.itertuples(index=False)]
    head = "| " + " | ".join(frame.columns) + " |\n|" + "---|" * len(frame.columns)
    return head + "\n" + "\n".join("| " + " | ".join(r) + " |" for r in cells)


def read(path: Path) -> np.ndarray:
    return np.asarray(tifffile.imread(path))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261105)
    checks, rows, morphs = [], [], []
    for crop in CROPS:
        gt = read(CROPS_DIR / crop / "labels.tif")
        for route, (ws_rel, config, label) in ROUTES.items():
            ws = RUNS / crop / ws_rel
            markers = read(ws / "normalized_markers" / f"{label}.tif")
            for flt in ("none", "500"):
                work = RUNS / crop / "phase1_closure" / f"{label}_{flt}"
                pred = read(work / "e09_dir" / "segmentations" / f"{config}.tif")
                ref = pd.read_csv(work / "evaluation" / "ignore_aware_results.csv").set_index("evaluation_protocol")
                base = {"crop_id": crop, "route": route, "size_filter": flt}

                a, _, morph = evaluate(pred, gt, markers, rule="dominant")
                # mean_matched_iou025 excluded: e09b's value depends on label order (see evaluate_v2 docstring)
                ok_a, d_a = same(a, ref.loc["seed_ignore"], skip=("mean_matched_iou025",))
                d_a025 = a["mean_matched_iou025"] - float(ref.loc["seed_ignore", "mean_matched_iou025"])
                perm = np.concatenate([[0], rng.permutation(np.arange(1, int(pred.max()) + 1)) + 1000])
                b, _, _ = evaluate(perm[pred.astype(np.int64)], gt, markers, rule="dominant", morph_iou=None)
                ok_b, d_b = same(b, pd.Series(a))
                c, _, _ = evaluate(pred, gt, None, rule="fraction", ignore_overlap_threshold=0.50, morph_iou=None)
                ok_c, d_c = same(c, ref.loc["overlap_ignore"], skip=("mean_matched_iou025",))
                d, _, _ = evaluate(pred, gt, markers, rule="fraction", ignore_overlap_threshold=0.50, morph_iou=None)
                ok_d, d_d = same(d, ref.loc["seed_ignore"], skip=("mean_matched_iou025",))
                checks.append(base | {
                    "A_markers_dominant_eq_e09b_seed_ignore": ok_a, "A_max_abs_diff": d_a,
                    "A_mean_matched_iou025_v2_minus_e09b": d_a025,
                    "B_permuted_labels_eq_A": ok_b, "B_max_abs_diff": d_b,
                    "C_no_markers_fraction050_eq_e09b_overlap_ignore": ok_c, "C_max_abs_diff": d_c,
                    "D_markers_fraction050_eq_e09b_seed_ignore": ok_d,
                    "e09b_seed_excluded_prediction_count": int(ref.loc["seed_ignore", "excluded_prediction_count"]),
                    "v2_seed_n_pred_ignored": a["n_pred_ignored"], "v2_seed_n_unseeded": a["n_pred_unseeded"],
                    "v2_seed_n_background_support": a["n_pred_background_support"],
                    "f1_050_e09b_seed": float(ref.loc["seed_ignore", "f1_iou050"]), "f1_050_v2_seed": a["f1_iou050"],
                    "f1_050_v2_no_markers": c["f1_iou050"], "f1_075_v2_no_markers": c["f1_iou075"]})
                for mode, r in (("markers_dominant", a), ("no_markers_fraction050", c), ("markers_fraction050", d)):
                    rows.append(base | {"mode": mode} | r)
                if not morph.empty:
                    morphs.append(morph.assign(**base))
                print(f"[v2] {crop} {route:18s} {flt:4s} A={ok_a} B={ok_b} C={ok_c} D={ok_d} "
                      f"F1@.50 e09b={ref.loc['seed_ignore', 'f1_iou050']:.6f} v2={a['f1_iou050']:.6f} "
                      f"nomk={c['f1_iou050']:.6f}", flush=True)

    checks = pd.DataFrame(checks)
    per = pd.DataFrame(rows)
    checks.to_csv(OUT / "verification.csv", index=False)
    per.to_csv(OUT / "per_crop_all_modes.csv", index=False)
    morph = pd.concat(morphs, ignore_index=True)
    morph.to_csv(OUT / "morphology_matched_seed.csv.gz", index=False)  # per object, ~7 MB: not tracked
    msum = morph.groupby(["route", "size_filter"], sort=False).agg(
        n_matched=("iou", "size"), volume_rel_error_median=("volume_rel_error", "median"),
        volume_rel_error_median_abs=("volume_rel_error", lambda x: x.abs().median()),
        d_eq_error_median=("d_eq_error", "median"), centroid_error_median=("centroid_error", "median"),
        centroid_error_p95=("centroid_error", lambda x: x.quantile(0.95)),
        surface_area_rel_error_median=("surface_area_rel_error", "median"),
        sphericity_error_median=("sphericity_error", "median"),
        elongation_error_median=("elongation_error", "median"),
        gt_d_eq_median=("gt_d_eq", "median"), gt_sphericity_median=("gt_sphericity", "median")).reset_index()
    msum.to_csv(OUT / "morphology_summary_pooled.csv", index=False)

    keys = ["route", "size_filter", "mode"]
    macro, pooled = [], []
    for key, g in per.groupby(keys, sort=False):
        m = dict(zip(keys, key), n_crops=len(g))
        p = dict(m)
        for s in ("025", "050", "075"):
            f1 = g[f"f1_iou{s}"]
            m |= {f"f1_iou{s}_mean": f1.mean(), f"f1_iou{s}_sd_ddof1": f1.std(ddof=1),
                  f"f1_iou{s}_median": f1.median(), f"f1_iou{s}_min": f1.min()}
            tp, fp, fn = (int(g[f"{x}_iou{s}"].sum()) for x in ("tp", "fp", "fn"))
            pr, rc = tp / (tp + fp), tp / (tp + fn)
            p |= {f"tp_iou{s}": tp, f"fp_iou{s}": fp, f"fn_iou{s}": fn, f"precision_iou{s}": pr,
                  f"recall_iou{s}": rc, f"f1_iou{s}": 2 * pr * rc / (pr + rc), f"ap_iou{s}": tp / (tp + fp + fn)}
        m |= {"ap_iou050_mean": g["ap_iou050"].mean(), "ap_iou075_mean": g["ap_iou075"].mean(),
              "mean_matched_iou050_mean": g["mean_matched_iou050"].mean()}
        macro.append(m)
        pooled.append(p)
    macro, pooled = pd.DataFrame(macro), pd.DataFrame(pooled)
    macro.to_csv(OUT / "crop_macro_by_mode.csv", index=False)
    pooled.to_csv(OUT / "pooled_object_by_mode.csv", index=False)

    n = len(checks)
    lines = ["# evaluate_v2 verification against frozen e09b (Phase 1d)", "",
             f"Runs: 12 crops x 3 routes x size filter (none, 500) = {n} evaluations.", "",
             "| check | passing |", "|---|---|"]
    for col in ("A_markers_dominant_eq_e09b_seed_ignore", "B_permuted_labels_eq_A",
                "C_no_markers_fraction050_eq_e09b_overlap_ignore", "D_markers_fraction050_eq_e09b_seed_ignore"):
        lines.append(f"| {col} | {int(checks[col].sum())}/{n} |")
    d25 = checks["A_mean_matched_iou025_v2_minus_e09b"]
    lines += ["", "A/C/D compare counts, P, R, F1 at IoU 0.25/0.50/0.75 and mean matched IoU at 0.50/0.75 (exact; floats",
              "within 1e-12). mean matched IoU at 0.25 is excluded: e09b keeps an arbitrary maximum-cardinality matching",
              "(label-order dependent); evaluate_v2 takes the max-total-IoU one among them. v2 - e09b at 0.25:",
              f"min {d25.min():+.6f}, max {d25.max():+.6f}, nonzero in {int((d25.abs() > 1e-12).sum())}/{n} runs.",
              "B compares everything, including mean matched IoU at 0.25."]
    lines += ["", "## Crop-macro F1 (sample SD), by mode", "",
              macro[keys + ["f1_iou050_mean", "f1_iou050_sd_ddof1", "f1_iou050_median", "f1_iou050_min",
                            "f1_iou075_mean", "f1_iou025_mean", "ap_iou050_mean"]].pipe(md),
              "", "## Pooled object, by mode", "",
              pooled[keys + ["tp_iou050", "fp_iou050", "fn_iou050", "f1_iou050", "f1_iou075", "ap_iou050"]].pipe(md),
              "", "## Morphology of objects matched at IoU 0.50 (markers_dominant; pooled objects, voxel units)", "",
              msum.pipe(md)]
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    hard = checks["A_markers_dominant_eq_e09b_seed_ignore"].all() and checks["B_permuted_labels_eq_A"].all()
    return 0 if hard else 1


if __name__ == "__main__":
    raise SystemExit(main())
