#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT = Path(".")

BASE_WITH_LS = ROOT / "results_batos/paper_tables/final_batos_ls/table_final_comparison_with_batos_ls.csv"
AUTO = ROOT / "results_batos/paper_tables/batos_local_split_auto/local_split_auto_summary.csv"

OUT = ROOT / "results_batos/paper_tables/final_batos_ls_auto"
OUT.mkdir(parents=True, exist_ok=True)

if not BASE_WITH_LS.exists():
    raise FileNotFoundError(f"Missing: {BASE_WITH_LS}")

if not AUTO.exists():
    raise FileNotFoundError(f"Missing: {AUTO}")

base = pd.read_csv(BASE_WITH_LS)
auto = pd.read_csv(AUTO)

a = auto.iloc[0]

auto_row = {
    "method": "BA-TOS-LS-AUTO",
    "mean_precision": a["mean_precision"],
    "mean_recall": a["mean_recall"],
    "F1@0.25": a["F1@0.25"],
    "F1@0.50": a["F1@0.50"],
    "F1@0.75": a["F1@0.75"],
    "mean_overseg": a["mean_overseg"],
    "total_fp": int(a["total_fp"]),
    "total_fn": int(a["total_fn"]),
}

df = pd.concat([base, pd.DataFrame([auto_row])], ignore_index=True)

order = {
    "Min-tree-core+WS": 0,
    "Otsu+WS": 1,
    "BA-TOS": 2,
    "BA-TOS-LS": 3,
    "BA-TOS-LS-AUTO": 4,
}

df["order"] = df["method"].map(order)
df = df.sort_values("order").drop(columns=["order"])

csv_out = OUT / "table_final_comparison_with_batos_ls_auto.csv"
tex_out = OUT / "table_final_comparison_with_batos_ls_auto_latex.tex"

df.to_csv(csv_out, index=False)


def fmt(x):
    if pd.isna(x):
        return "--"
    return f"{float(x):.3f}"


latex = []
latex.append(r"\begin{table}[t]")
latex.append(r"\centering")
latex.append(r"\caption{Instance-level comparison including manual and automatic local-split refinements.}")
latex.append(r"\label{tab:main_results_batos_ls_auto}")
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
