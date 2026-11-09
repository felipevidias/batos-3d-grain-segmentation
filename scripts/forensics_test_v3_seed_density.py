#!/usr/bin/env python3
"""Phase 5e step 4, supplementary forensics (diagnosis only; nothing tuned, frozen v3 unchanged).

The frozen aggregator's first-stage rule (scripts/aggregate_test_v3.py) compares ToS markers with the VALID GT count
only. ToS markers seed every grain, including the face-touching ones excluded from scoring, so the relevant density is
markers / (valid + excluded GT). This script adds that density and the candidate-capacity saturation per crop.

Inputs : results/test_v3/{per_crop_all_modes,forensics_stages,crop_info}.csv (from aggregate_test_v3.py)
Outputs: results/test_v3/seed_density.csv (per crop), seed_density_by_class.csv, seed_density.md
GT enters only as counts (diagnosis), never as a parameter.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "test_v3"
K_MIN_V2 = 1800  # MIN half of true_dual --max-selected-nodes 3600; v3 rule v2/s^3 (src/batos_v3/constants.py)


def main() -> int:
    per = pd.read_csv(OUT / "per_crop_all_modes.csv", dtype={"size_filter": str})
    x = per[(per.method == "v3") & (per.route == "BA-TOS_sel") & (per.size_filter == "none")
            & (per["mode"] == "markers_dominant")].set_index("crop_id")
    f = pd.read_csv(OUT / "forensics_stages.csv").set_index("crop_id")
    info = pd.read_csv(OUT / "crop_info.csv").set_index("crop_id")
    d = pd.DataFrame({
        "specimen": f.specimen, "test_class": f.test_class, "s": info.v3_s,
        "k_min_cap": (K_MIN_V2 / info.v3_s ** 3).round().astype(int), "min_candidates": f.min_candidates,
        "min_in_internal": f.min_in_internal, "tos_markers": f.tos_markers,
        "gt_valid": x.valid_gt_count, "gt_excluded": x.n_gt_excluded,
        "pred_total": x.n_pred_total, "pred_ignored": x.n_pred_ignored, "pred_evaluable": x.prediction_count_evaluable,
        "v3_f1_050": x.f1_iou050, "fixflat_f1_050": f.fixflat_f1_sel, "flag": f.flag})
    d["gt_all"] = d.gt_valid + d.gt_excluded
    d["cap_saturated"] = d.min_candidates >= d.k_min_cap
    d["candidates_per_grain"] = d.min_candidates / d.gt_all
    d["markers_per_grain"] = d.tos_markers / d.gt_all
    d.to_csv(OUT / "seed_density.csv")

    g = d.groupby(["specimen", "test_class"]).agg(
        n=("v3_f1_050", "size"), s=("s", "median"), k_min_cap=("k_min_cap", "median"),
        saturated=("cap_saturated", "sum"), gt_valid=("gt_valid", "median"), gt_all=("gt_all", "median"),
        candidates_per_grain=("candidates_per_grain", "median"), min_in_internal=("min_in_internal", "median"),
        markers_per_grain=("markers_per_grain", "median"), v3_f1_median=("v3_f1_050", "median"),
        n_flag=("flag", "sum")).reset_index()
    g.to_csv(OUT / "seed_density_by_class.csv", index=False)

    L = ["| specimen | class | n | s | K_MIN cap (1800/s³) | cap saturated | GT valid / all (median) "
         "| candidates / grain | MIN in internal | ToS markers / grain | v3 F1@0.50 median | flagged |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in g.itertuples():
        L.append(f"| {r.specimen} | {r.test_class} | {r.n} | {r.s:g} | {r.k_min_cap:.0f} | {r.saturated}/{r.n} "
                 f"| {r.gt_valid:.0f} / {r.gt_all:.0f} | {r.candidates_per_grain:.2f} | {r.min_in_internal:.0f} "
                 f"| {r.markers_per_grain:.2f} | {r.v3_f1_median:.3f} | {r.n_flag} |")
    L += ["", "| crop | class | candidates / grain | ToS markers / grain | pred evaluable / GT valid | v3 F1@0.50 | fixflat |",
          "|---|---|---|---|---|---|---|"]
    for c, r in d[d.flag].iterrows():
        L.append(f"| {c} | {r.test_class} | {r.candidates_per_grain:.2f} | {r.markers_per_grain:.2f} "
                 f"| {r.pred_evaluable} / {r.gt_valid} | {r.v3_f1_050:.3f} | {r.fixflat_f1_050:.3f} |")
    (OUT / "seed_density.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
