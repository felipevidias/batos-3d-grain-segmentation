#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

OUT = Path("results_batos/paper_tables/final_revision")
OUT.mkdir(parents=True, exist_ok=True)

rows = []

# EFRGP01_00: BA-TOS + Otsu
p00 = Path("results_batos/paper_tables/batos_min500_vs_otsu_selected_crops.csv")
if p00.exists():
    df = pd.read_csv(p00)
    for _, r in df.iterrows():
        rows.append({
            "scan": "EFRGP01_00",
            "crop_id": r["crop_id"],
            "method": r["method"],
            "precision050": float(r["precision050"]),
            "recall050": float(r["recall050"]),
            "f1050": float(r["f1050"]),
            "fp050": int(r["fp050"]),
            "fn050": int(r["fn050"]),
            "overseg_ratio": float(r["overseg_ratio"]),
        })
else:
    print("MISSING:", p00)

# EFRGP01_00: Min-tree-core+WS
p00m = Path("results_batos/baselines/macke_mintree_selected_crops/macke_mintree_selected_crops_summary.csv")
if p00m.exists():
    df = pd.read_csv(p00m)
    for _, r in df.iterrows():
        rows.append({
            "scan": "EFRGP01_00",
            "crop_id": r["crop_id"],
            "method": "Min-tree-core+WS",
            "precision050": float(r["precision050"]),
            "recall050": float(r["recall050"]),
            "f1050": float(r["f1050"]),
            "fp050": int(r["fp050"]),
            "fn050": int(r["fn050"]),
            "overseg_ratio": float(r["overseg_ratio"]),
        })
else:
    print("MISSING:", p00m)

# EFRGP01_01: all methods
p01 = Path("results_batos/EFRGP01_01_selected_crops/EFRGP01_01_selected_crops_methods_summary.csv")
if p01.exists():
    df = pd.read_csv(p01)
    for _, r in df.iterrows():
        method = r["method"]
        if method == "Macke-MinTree+WS":
            method = "Min-tree-core+WS"

        rows.append({
            "scan": "EFRGP01_01",
            "crop_id": r["crop_id"],
            "method": method,
            "precision050": float(r["precision050"]),
            "recall050": float(r["recall050"]),
            "f1050": float(r["f1050"]),
            "fp050": int(r["fp050"]),
            "fn050": int(r["fn050"]),
            "overseg_ratio": float(r["overseg_ratio"]),
        })
else:
    print("MISSING:", p01)

all_df = pd.DataFrame(rows)

methods = ["Min-tree-core+WS", "Otsu+WS", "BA-TOS"]
all_df = all_df[all_df["method"].isin(methods)].copy()

# sanity check
print("\nCounts by method:")
print(all_df.groupby("method").size())

all_csv = OUT / "primary_existing_summaries_all_crops.csv"
all_df.to_csv(all_csv, index=False)

summary = (
    all_df
    .groupby("method", as_index=False)
    .agg(
        mean_precision=("precision050", "mean"),
        mean_recall=("recall050", "mean"),
        mean_f1=("f1050", "mean"),
        mean_overseg=("overseg_ratio", "mean"),
        total_fp=("fp050", "sum"),
        total_fn=("fn050", "sum"),
        n=("f1050", "count"),
    )
)

order = {m:i for i,m in enumerate(methods)}
summary["order"] = summary["method"].map(order)
summary = summary.sort_values("order").drop(columns=["order"])

summary_csv = OUT / "table_primary_results_iou050_FROM_EXISTING_SUMMARIES.csv"
summary.to_csv(summary_csv, index=False)

def fmt(x):
    return f"{float(x):.3f}"

latex = []
latex.append(r"\begin{table}[t]")
latex.append(r"\centering")
latex.append(r"\caption{Instance-level comparison across two consecutive scans and twelve representative 3D crops at IoU@0.50.}")
latex.append(r"\label{tab:main_results_iou050}")
latex.append(r"\resizebox{\columnwidth}{!}{%")
latex.append(r"\begin{tabular}{lrrrrrr}")
latex.append(r"\hline")
latex.append(r"Method & P@0.50 & R@0.50 & F1@0.50 & Overseg. & FP & FN \\")
latex.append(r"\hline")
for _, r in summary.iterrows():
    latex.append(
        f"{r['method']} & {fmt(r['mean_precision'])} & {fmt(r['mean_recall'])} & {fmt(r['mean_f1'])} & {fmt(r['mean_overseg'])} & {int(r['total_fp'])} & {int(r['total_fn'])} \\\\"
    )
latex.append(r"\hline")
latex.append(r"\end{tabular}%")
latex.append(r"}")
latex.append(r"\end{table}")

latex_path = OUT / "table_primary_results_iou050_FROM_EXISTING_SUMMARIES_latex.tex"
latex_path.write_text("\n".join(latex), encoding="utf-8")

print("\nwrote:", all_csv)
print("wrote:", summary_csv)
print("wrote:", latex_path)
print("\n=== FINAL PRIMARY TABLE ===")
print(summary.to_string(index=False))
