#!/usr/bin/env python3
"""Phase 5a: rank parameters by effect size from results/sensitivity/global_sensitivity.csv.

Effect size of a parameter = max over its perturbations of |delta crop-macro mean F1@0.50| (12 crops,
evaluate_v2 with markers == e09b seed_ignore), per route. Influential = effect > 0.01 on either perturbation
(primary route BA-TOS_sel; ToS reported alongside). A perturbation whose pipeline fails on any crop is
reported as `breaks` (not as an effect size). Evaluator rows (kind eval) are protocol, not method, and are
excluded from the parameter count. Writes results/sensitivity/{ranking.csv, RANKING.md}.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "results" / "sensitivity" / "global_sensitivity.csv"
OUT = REPO / "results" / "sensitivity"
THRESHOLD = 0.01


def md(frame: pd.DataFrame) -> str:
    cells = [[f"{v:+.4f}" if isinstance(v, float) and c.startswith("delta") else
              (f"{v:.4f}" if isinstance(v, float) else str(v)) for c, v in zip(frame.columns, r)]
             for r in frame.itertuples(index=False)]
    head = "| " + " | ".join(frame.columns) + " |\n|" + "---|" * len(frame.columns)
    return head + "\n" + "\n".join("| " + " | ".join(r) + " |" for r in cells)


def main() -> int:
    g = pd.read_csv(SRC)
    ident = g[g["param"] == "identity"]
    g = g[g["param"] != "identity"]
    rows = []
    for (group, param), sub in g.groupby(["group", "param"], sort=False):
        row = {"group": group, "param": param, "kind": sub["kind"].iloc[0], "stage": sub["stage"].iloc[0],
               "frozen_value": sub["frozen_value"].iloc[0]}
        complete = True
        for route, tag in (("BA-TOS_sel", "sel"), ("BA-TOS_ToS", "tos")):
            r = sub[sub["route"] == route]
            failed = r["failed_crops"].fillna("").astype(str).str.len() > 0
            incomplete = r["mean_f1_050"].isna() & ~failed
            complete &= not incomplete.any()
            d = r["delta_mean_f1_050"].where(~failed)
            dm = r["delta_min_f1_050"].where(~failed)
            pert = ", ".join(f"{a}={b}: {c:+.4f}" if np.isfinite(c) else f"{a}={b}: {'FAIL' if f else 'n/a'}"
                             for a, b, c, f in zip(r["direction"], r["perturbed_value"], d, failed))
            row |= {f"{tag}_max_abs_delta_mean": float(d.abs().max()) if d.notna().any() else np.nan,
                    f"{tag}_max_abs_delta_min": float(dm.abs().max()) if dm.notna().any() else np.nan,
                    f"{tag}_perturbations": pert, f"{tag}_breaks": bool(failed.any())}
        row["complete"] = complete
        row["influential_sel"] = bool(row["sel_max_abs_delta_mean"] > THRESHOLD)
        row["influential_any_route"] = bool(row["influential_sel"] or row["tos_max_abs_delta_mean"] > THRESHOLD)
        rows.append(row)
    rank = pd.DataFrame(rows).sort_values(["sel_max_abs_delta_mean", "tos_max_abs_delta_mean"], ascending=False,
                                          na_position="last")
    rank.to_csv(OUT / "ranking.csv", index=False)
    method = rank[rank["kind"] != "eval"]
    method_done = method[method["complete"]]
    infl = method_done[method_done["influential_any_route"]]
    lines = ["# Phase 5a — parameter ranking by effect size", "",
             f"Source: `global_sensitivity.csv`. Effect = max |Δ crop-macro mean F1@0.50| over the perturbations "
             f"(12 crops, seed-ignore). Influential: > {THRESHOLD} on either perturbation, either route "
             "(BA-TOS_sel primary).", ""]
    if len(ident):
        lines += ["Identity run (whole pipeline re-run from scratch by the sensitivity driver):", "",
                  md(ident[["route", "n_crops_done", "mean_f1_050", "min_f1_050", "delta_mean_f1_050",
                            "delta_min_f1_050"]]), ""]
    lines += [f"Method parameters swept and complete: {len(method_done)} / {len(method)}; "
              f"influential: **{len(infl)}** "
              f"(sel route only: {int(method_done['influential_sel'].sum())}).", "",
              "## Influential parameters", "",
              md(infl[["group", "param", "frozen_value", "sel_max_abs_delta_mean", "tos_max_abs_delta_mean",
                       "sel_perturbations"]]) if len(infl) else "(none yet)", "",
              "## All parameters (ranked)", "",
              md(rank[["group", "param", "kind", "frozen_value", "sel_max_abs_delta_mean", "sel_max_abs_delta_min",
                       "tos_max_abs_delta_mean", "sel_breaks", "complete"]])]
    (OUT / "RANKING.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[:12]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
