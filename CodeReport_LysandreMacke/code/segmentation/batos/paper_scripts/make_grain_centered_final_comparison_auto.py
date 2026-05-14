#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT = Path(".")

ORIG = ROOT / "results_batos/paper_tables/grain_centered_analysis/grain_centered_summary_by_method.csv"
LS = ROOT / "results_batos/paper_tables/grain_centered_analysis_batos_ls/grain_centered_summary_by_method.csv"
AUTO = ROOT / "results_batos/paper_tables/grain_centered_analysis_batos_ls_auto/grain_centered_summary_by_method.csv"

OUT = ROOT / "results_batos/paper_tables/grain_centered_final_auto"
OUT.mkdir(parents=True, exist_ok=True)

paths = [ORIG, LS, AUTO]

for p in paths:
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")

dfs = []

for p in paths:
    df = pd.read_csv(p)
    dfs.append(df)

final = pd.concat(dfs, ignore_index=True)

order = {
    "Min-tree-core+WS": 0,
    "Otsu+WS": 1,
    "BA-TOS": 2,
    "BA-TOS-LS": 3,
    "BA-TOS-LS-AUTO": 4,
}

final["order"] = final["method"].map(order)
final = final.sort_values("order").drop(columns=["order"])

csv_out = OUT / "grain_centered_final_comparison_auto.csv"
tex_out = OUT / "grain_centered_final_comparison_auto_latex.tex"

final.to_csv(csv_out, index=False)


def fmt(x):
    return f"{float(x):.3f}"


latex = []
latex.append(r"\begin{table}[t]")
latex.append(r"\centering")
latex.append(r"\caption{Grain-centered evaluation over valid reference grains.}")
latex.append(r"\label{tab:grain_centered_auto}")
latex.append(r"\resizebox{\columnwidth}{!}{%")
latex.append(r"\begin{tabular}{lrrrrrr}")
latex.append(r"\hline")
latex.append(r"Method & Rec@.50 & Rec@.75 & Mean IoU & Median IoU & Med. Centroid Err. & Med. Vol. Err. \\")
latex.append(r"\hline")

for _, r in final.iterrows():
    latex.append(
        f"{r['method']} & "
        f"{fmt(r['recall_iou050'])} & "
        f"{fmt(r['recall_iou075'])} & "
        f"{fmt(r['mean_best_iou'])} & "
        f"{fmt(r['median_best_iou'])} & "
        f"{fmt(r['median_centroid_error'])} & "
        f"{fmt(r['median_relative_volume_error'])} \\\\"
    )

latex.append(r"\hline")
latex.append(r"\end{tabular}%")
latex.append(r"}")
latex.append(r"\end{table}")

tex_out.write_text("\n".join(latex), encoding="utf-8")

print(final.to_string(index=False))
print()
print("wrote:", csv_out)
print("wrote:", tex_out)
