#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT = Path(".")

BASE = ROOT / "results_batos/paper_tables/final_revision_clean/table_manual_clean_combined_results.csv"
LS_MULTI = ROOT / "results_batos/paper_tables/batos_local_split_large_labels/local_split_multi_iou_summary.csv"

OUT = ROOT / "results_batos/paper_tables/final_batos_ls"
OUT.mkdir(parents=True, exist_ok=True)

if not BASE.exists():
    raise FileNotFoundError(f"Arquivo base não encontrado: {BASE}")

if not LS_MULTI.exists():
    raise FileNotFoundError(
        f"Resumo multi-IoU do BA-TOS-LS não encontrado: {LS_MULTI}\n"
        "Rode antes: python code/segmentation/batos/paper_scripts/evaluate_batos_ls_multi_iou_from_manifest.py"
    )

base = pd.read_csv(BASE)
ls_multi = pd.read_csv(LS_MULTI)

ls_row = ls_multi.iloc[0]

local = {
    "method": "BA-TOS-LS",
    "mean_precision": ls_row["mean_precision"],
    "mean_recall": ls_row["mean_recall"],
    "F1@0.25": ls_row["F1@0.25"],
    "F1@0.50": ls_row["F1@0.50"],
    "F1@0.75": ls_row["F1@0.75"],
    "mean_overseg": ls_row["mean_overseg"],
    "total_fp": int(ls_row["total_fp"]),
    "total_fn": int(ls_row["total_fn"]),
}

rows = []

for _, r in base.iterrows():
    rows.append({
        "method": r["method"],
        "mean_precision": r["mean_precision"],
        "mean_recall": r["mean_recall"],
        "F1@0.25": r["F1@0.25"],
        "F1@0.50": r["F1@0.50"],
        "F1@0.75": r["F1@0.75"],
        "mean_overseg": r["mean_overseg"],
        "total_fp": int(r["total_fp"]),
        "total_fn": int(r["total_fn"]),
    })

rows.append(local)

df = pd.DataFrame(rows)

order = {
    "Min-tree-core+WS": 0,
    "Otsu+WS": 1,
    "BA-TOS": 2,
    "BA-TOS-LS": 3,
}

df["order"] = df["method"].map(order)
df = df.sort_values("order").drop(columns=["order"])

csv_out = OUT / "table_final_comparison_with_batos_ls.csv"
tex_out = OUT / "table_final_comparison_with_batos_ls_latex.tex"

df.to_csv(csv_out, index=False)


def fmt(x):
    if pd.isna(x):
        return "--"
    return f"{float(x):.3f}"


latex = []
latex.append(r"\begin{table}[t]")
latex.append(r"\centering")
latex.append(r"\caption{Instance-level comparison including the local-split refinement.}")
latex.append(r"\label{tab:main_results_batos_ls}")
latex.append(r"\resizebox{\columnwidth}{!}{%")
latex.append(r"\begin{tabular}{lrrrrrrrr}")
latex.append(r"\hline")
latex.append(r"Method & P@.50 & R@.50 & F1@.25 & F1@.50 & F1@.75 & Overseg. & FP & FN \\")
latex.append(r"\hline")

for _, r in df.iterrows():
    latex.append(
        f"{r['method']} & "
        f"{fmt(r['mean_precision'])} & "
        f"{fmt(r['mean_recall'])} & "
        f"{fmt(r['F1@0.25'])} & "
        f"{fmt(r['F1@0.50'])} & "
        f"{fmt(r['F1@0.75'])} & "
        f"{fmt(r['mean_overseg'])} & "
        f"{int(r['total_fp'])} & "
        f"{int(r['total_fn'])} \\\\"
    )

latex.append(r"\hline")
latex.append(r"\end{tabular}%")
latex.append(r"}")
latex.append(r"\end{table}")

tex_out.write_text("\n".join(latex), encoding="utf-8")

print(df.to_string(index=False))
print()
print("wrote:", csv_out)
print("wrote:", tex_out)
