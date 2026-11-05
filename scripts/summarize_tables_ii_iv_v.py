#!/usr/bin/env python3
"""Phase 1e: compare the frozen-code reproduction with paper Tables II, IV and V.

Reads C:\\BA_TOS_DATA\\runs\\reproduction_tables\\{mintree,sensitivity_crop04,crop_size} (batos_v2.reproduce_tables)
and results/batos_sibgrapi/crop_macro_summary.csv (Phase 1 closure). Writes to results/batos_sibgrapi/:
  table2_mintree_per_crop.csv, table2_full_comparison.csv, table4_sensitivity_comparison.csv,
  table5_crop_size_comparison.csv, tables_II_IV_V_reproduction.md
and appends phase1e rows to results/all_runs.csv (rows with experiment_id starting "phase1e/" are replaced).
Primary protocol: seed_ignore. Metrics are never adjusted; every difference from the paper is listed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
ROOT = Path(r"C:\BA_TOS_DATA\runs\reproduction_tables")
FROZEN = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen")
OUT = REPO / "results" / "batos_sibgrapi"
CROPS = [f"EFRGP01_0{s}_crop0{c}" for s in (0, 1) for c in (2, 3, 4, 5, 6, 7)]
TOL = 0.0005  # rounding tolerance of 3-decimal paper values

PAPER_T2 = {
    "BA-TOS_ToS": {"mean": 0.824, "sd": 0.367, "median": 0.985, "min": 0.012, "mean_075": 0.817},
    "BA-TOS_sel": {"mean": 0.906, "sd": 0.210, "median": 0.985, "min": 0.260, "mean_075": 0.899},
    "Otsu+WS(edt_otsu)": {"mean": 0.830, "sd": 0.182, "median": 0.882, "min": 0.260, "mean_075": 0.770},
    "Min-tree-core+WS(mintree_e14)": {"mean": 0.768, "sd": 0.315, "median": 0.890, "min": 0.088, "mean_075": 0.751},
}
PAPER_T4 = {("score_threshold", v): 0.985 for v in (0.25, 0.30, 0.35, 0.40, 0.45)} | {
    ("k", 900): 0.510, ("k", 1800): 0.833, ("k", 3600): 0.985, ("k", 7200): 0.998,
    ("nms_radius", 4): 0.976, ("nms_radius", 6): 0.985, ("nms_radius", 8): 0.985, ("nms_radius", 10): 0.985}
PAPER_T5 = {("EFRGP01_00_crop04", 128): 0.830, ("EFRGP01_00_crop04", 160): 0.979, ("EFRGP01_00_crop04", 200): 0.985,
            ("EFRGP01_01_crop04", 128): 0.871, ("EFRGP01_01_crop04", 160): 0.955, ("EFRGP01_01_crop04", 200): 0.985}


def verdict(diff: float) -> str:
    return "match" if abs(diff) <= TOL else "DIFFERS"


def load(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def seed(results: list[dict]) -> pd.DataFrame:
    frame = pd.DataFrame(results)
    return frame[frame.evaluation_protocol == "seed_ignore"]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    runs_rows = []

    # ---- Table II: Min-tree per crop + consistency of the ToS marker re-run in the same watershed call
    per, consistency = [], []
    for crop in CROPS:
        rec = load(ROOT / "mintree" / crop / "run.json")
        s = seed(rec["e09b_results"]).set_index("marker_source")
        frozen = seed(load(FROZEN / crop / "frozen_run.json")["e09b_results"]).set_index("marker_source")
        m = s.loc["mintree_e14"]
        per.append({"crop_id": crop, **{k: m[k] for k in (
            "valid_gt_count", "prediction_count_evaluable", "tp_iou050", "fp_iou050", "fn_iou050",
            "precision_iou050", "recall_iou050", "f1_iou025", "f1_iou050", "f1_iou075")}})
        consistency.append(abs(float(s.loc["tos_e11_min", "f1_iou050"]) - float(frozen.loc["tos_e11_min", "f1_iou050"])))
        runs_rows.append(run_row(f"phase1e/table2_mintree/{crop}", "Min-tree-core+WS(mintree_e14)", crop, rec, m,
                                 "E14 mintree K3600 area 20-45000; NMS=6; size_filter=none; protocol=seed_ignore; "
                                 "runtime excludes shared preprocessing; concurrent_jobs=3"))
    per = pd.DataFrame(per)
    per.to_csv(OUT / "table2_mintree_per_crop.csv", index=False)
    if max(consistency) != 0.0:
        raise RuntimeError(f"tos_e11_min in the Min-tree watershed differs from the frozen run: {consistency}")

    macro = pd.read_csv(OUT / "crop_macro_summary.csv")
    stats = {}
    for route in ("BA-TOS_ToS", "BA-TOS_sel", "Otsu+WS(edt_otsu)"):
        r = macro[(macro.route == route) & (macro.size_filter == "none")].iloc[0]
        stats[route] = {"mean": r.f1_050_mean, "sd": r.f1_050_sd_sample, "median": r.f1_050_median,
                        "min": r.f1_050_min, "mean_075": r.f1_075_mean}
    f = per.f1_iou050.astype(float)
    stats["Min-tree-core+WS(mintree_e14)"] = {"mean": f.mean(), "sd": f.std(ddof=1), "median": f.median(),
                                              "min": f.min(), "mean_075": per.f1_iou075.astype(float).mean()}
    t2 = pd.DataFrame([{"route": route, "statistic": k, "paper": PAPER_T2[route][k], "ours": round(float(v), 4),
                        "diff": round(float(v) - PAPER_T2[route][k], 4)}
                       for route, d in stats.items() for k, v in d.items()])
    t2["verdict"] = t2["diff"].map(verdict)
    t2.to_csv(OUT / "table2_full_comparison.csv", index=False)

    # ---- Table IV
    rec = load(ROOT / "sensitivity_crop04" / "run.json")
    sens = pd.read_csv(ROOT / "sensitivity_crop04" / "development_sensitivity_results.csv")
    t4 = []
    for _, r in sens.iterrows():
        for factor in ("score_threshold", "k", "nms_radius"):
            default = {"score_threshold": 0.35, "k": 3600, "nms_radius": 6}
            if all(r[o] == v for o, v in default.items() if o != factor):
                key = (factor, float(r[factor]) if factor == "score_threshold" else int(r[factor]))
                if key in PAPER_T4:
                    t4.append({"factor": factor, "value": key[1], "label": r.label, "paper_f1_050": PAPER_T4[key],
                               "ours_f1_050": round(float(r.f1_iou050), 4),
                               "diff": round(float(r.f1_iou050) - PAPER_T4[key], 4),
                               "ours_f1_075": round(float(r.f1_iou075), 4),
                               "n_markers_pred_eval": r.prediction_count_evaluable})
        steps = [s for s in rec["steps"] if s["step"].startswith(r.label + ":")]
        runs_rows.append(run_row(f"phase1e/table4_sensitivity/{r.label}", "BA-TOS_ToS", "EFRGP01_00_crop04",
                                 {"seconds_total": sum(s["seconds"] for s in steps),
                                  "peak_memory_mb": max((s["peak_memory_mb"] for s in steps), default=None),
                                  "git": rec["git"], "_path": str(ROOT / "sensitivity_crop04" / "run.json")}, r,
                                 f"tau_s={r.score_threshold}; K={r.k}; NMS={r.nms_radius}; size_filter=none; "
                                 "protocol=seed_ignore; runtime = this variant's selector/spatial/marker steps "
                                 "(shared watershed+eval excluded); concurrent_jobs=3"))
    t4 = pd.DataFrame(t4).sort_values(["factor", "value"])
    t4["verdict"] = t4["diff"].map(verdict)
    t4.to_csv(OUT / "table4_sensitivity_comparison.csv", index=False)

    # ---- Table V
    cs = pd.read_csv(ROOT / "crop_size" / "crop_size_results.csv")
    cs = cs[(cs.evaluation_protocol == "seed_ignore") & (cs.marker_source == "tos_e11_min")]
    t5 = []
    for _, r in cs.iterrows():
        key = (r.parent_crop_id, int(r.crop_size))
        sub = load(ROOT / "crop_size" / "e11" / r.crop_id / "frozen_run.json")
        gate = sub.get("route", {}).get("gate", sub.get("route", {}))
        t5.append({"parent": key[0], "size": key[1], "paper_f1_050": PAPER_T5[key],
                   "ours_f1_050": round(float(r.f1_iou050), 4), "diff": round(float(r.f1_iou050) - PAPER_T5[key], 4),
                   "ours_f1_075": round(float(r.f1_iou075), 4), "valid_gt": r.valid_gt_count,
                   "gate_triggered": gate.get("gate_triggered"), "markers": gate.get("selected_tos_marker_rows")})
        runs_rows.append(run_row(f"phase1e/table5_crop_size/{r.crop_id}", "BA-TOS_ToS", r.crop_id,
                                 sub | {"_path": str(ROOT / "crop_size" / "e11" / r.crop_id / "frozen_run.json")}, r,
                                 f"centre crop {key[1]}^3 of {key[0]}; parent scale params frozen; size_filter=none; "
                                 "protocol=seed_ignore; concurrent_jobs=3"))
    t5 = pd.DataFrame(t5).sort_values(["parent", "size"])
    t5["verdict"] = t5["diff"].map(verdict)
    means = t5.groupby("size")[["paper_f1_050", "ours_f1_050"]].mean().round(4)
    t5.to_csv(OUT / "table5_crop_size_comparison.csv", index=False)

    # ---- all_runs.csv
    all_runs = REPO / "results" / "all_runs.csv"
    existing = pd.read_csv(all_runs)
    existing = existing[~existing.experiment_id.astype(str).str.startswith("phase1e/")]
    pd.concat([existing, pd.DataFrame(runs_rows)], ignore_index=True).to_csv(all_runs, index=False)

    write_report(t2, t4, t5, means, per)
    for name, frame in (("Table II", t2), ("Table IV", t4), ("Table V", t5)):
        print(name)
        print(frame.to_string(index=False))
    return 0


def run_row(experiment_id: str, method: str, crop: str, rec: dict, m, parameters: str) -> dict:
    return {
        "experiment_id": experiment_id, "method": method, "model": "batos-paper-v2-frozen (E-series verbatim)",
        "volume": crop[:10], "crop": crop, "device": "cpu (AMD Ryzen 7 5700X)", "parameters": parameters,
        "runtime_seconds": rec.get("seconds_total"), "peak_memory_mb": rec.get("peak_memory_mb"),
        "tp": m["tp_iou050"], "fp": m["fp_iou050"], "fn": m["fn_iou050"], "precision": m["precision_iou050"],
        "recall": m["recall_iou050"], "f1_025": m["f1_iou025"], "f1_050": m["f1_iou050"], "f1_075": m["f1_iou075"],
        "ap50": None, "ap75": None, "mean_matched_iou": None, "n_pred": m["prediction_count_evaluable"],
        "n_gt": m["valid_gt_count"], "git_commit": rec["git"]["commit"],
        "run_json": rec.get("_path", ""), "memory_run_json": rec.get("_path", ""),
    }


def md_table(frame: pd.DataFrame) -> str:
    cols = list(frame.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in frame.itertuples(index=False)]
    return "\n".join(lines)


def write_report(t2, t4, t5, means, per) -> None:
    diffs = pd.concat([t2.assign(table="II"), t4.assign(table="IV"), t5.assign(table="V")])
    n_diff = int((diffs.verdict == "DIFFERS").sum())
    text = f"""# Reproduction of paper Tables II, IV and V with the frozen E-series code

Generated by `scripts/summarize_tables_ii_iv_v.py` from runs of `python -I -m batos_v2.reproduce_tables`
(outputs in `C:\\BA_TOS_DATA\\runs\\reproduction_tables\\`) and the Phase 1 closure (`batos_v2.frozen`, 12 crops).
Protocol: seed_ignore, no prediction size filter, crop-macro SD = sample SD. "match" = |diff| <= {TOL}
(3-decimal rounding of the paper). **{n_diff} statistic(s) differ beyond rounding.**

## Table II (crop-macro F1, 12 crops)

{md_table(t2)}

Min-tree-core+WS = E14 Min-tree source control (`e14_build_mintree_candidates.py --max-candidates 3600 --area-min 20
--area-max 45000`, then the same E11 marker builder and fixed-mask watershed as BA-TOS). The ToS marker set was re-run
in the same watershed call on every crop and equals the frozen run exactly (consistency check). On EFRGP01_00_crop05
the Min-tree output equals the recovered original E13 output (all e09b metrics; scores within 1.7e-16).

Min-tree per crop (F1@0.50): {", ".join(f"{c[-11:]} {v:.3f}" for c, v in zip(per.crop_id, per.f1_iou050))}

## Table IV (sensitivity, EFRGP01_00_crop04, F1@0.50)

{md_table(t4)}

Flags of `run_e13_development_sensitivity.sh`: tau_s and K variants re-run the true-dual selector, spatialization and
marker builder; NMS variants reuse the default (tau 0.35, K 3600) spatial diagnostics. One fixed-mask watershed for all
variants, mask = default K3600/NMS6 markers + baseline markers.

## Table V (crop extent, centre crops of the two crop04 locations, BA-TOS_ToS F1@0.50)

{md_table(t5)}

Mean over the two locations:

{md_table(means.reset_index())}

`e15_make_crop_size_suite.py` (unmodified) cuts centred 128/160/200 cubes from each 200-cube parent and copies the parent
scale parameters (no re-estimation); each sub-crop then runs the full frozen pipeline (`batos_v2.frozen`).

## Notes
- Runtimes in `results/all_runs.csv` for these rows were measured with three jobs running concurrently (Min-tree,
  sensitivity, crop size) on 8 cores / 32 GB; they are upper bounds on the single-job runtime.
"""
    (OUT / "tables_II_IV_V_reproduction.md").write_text(text, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
