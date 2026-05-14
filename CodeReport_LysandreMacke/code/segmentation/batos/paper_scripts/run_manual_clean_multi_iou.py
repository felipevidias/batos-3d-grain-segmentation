#!/usr/bin/env python3
from pathlib import Path
import csv
import subprocess
import sys
import pandas as pd
import numpy as np

MANIFEST = Path("results_batos/paper_tables/final_revision_clean/manual_clean_manifest_multi_iou.csv")
OUT = Path("results_batos/paper_tables/final_revision_clean")
RUNS = OUT / "manual_metrics_runs"
INSTANCE = Path("code/segmentation/batos/batos_instance_metrics.py")

OUT.mkdir(parents=True, exist_ok=True)
RUNS.mkdir(parents=True, exist_ok=True)

if not MANIFEST.exists():
    raise SystemExit(f"missing manifest: {MANIFEST}")

if not INSTANCE.exists():
    raise SystemExit(f"missing instance metrics script: {INSTANCE}")

thresholds = ["0.25", "0.50", "0.75"]

manifest = list(csv.DictReader(open(MANIFEST)))
valid = [r for r in manifest if r["gt_found"] == "1" and r["pred_found"] == "1"]

print("valid rows:", len(valid), "/", len(manifest))

if len(valid) != len(manifest):
    raise SystemExit("Manifest has missing files. Fix before running metrics.")

rows = []

for r in valid:
    scan = r["scan"]
    crop = r["crop_id"]
    method = r["method"]

    safe_method = (
        method.replace("+", "plus")
        .replace("/", "_")
        .replace(" ", "_")
        .replace("-", "_")
    )

    for th in thresholds:
        safe_th = th.replace(".", "p")
        outdir = RUNS / scan / crop / safe_method / f"iou_{safe_th}"
        outdir.mkdir(parents=True, exist_ok=True)

        name = f"{scan}_{crop}_{safe_method}_iou{safe_th}"

        cmd = [
            sys.executable, str(INSTANCE),
            "--gt", r["gt"],
            "--pred", r["pred"],
            "--out-dir", str(outdir),
            "--min-gt-area", "1000",
            "--ignore-gt-border",
            "--thresholds", th,
            "--name", name,
        ]

        print("RUN:", " ".join(cmd))
        subprocess.run(cmd, check=True)

        metric_file = outdir / f"{name}_instance_metrics.csv"
        df = pd.read_csv(metric_file)

        for _, m in df.iterrows():
            rows.append({
                "scan": scan,
                "crop_id": crop,
                "method": method,
                "threshold": float(m["threshold"]),
                "gt_valid": int(m["gt_valid"]),
                "pred_total": int(m["pred_total"]),
                "pred_eval": int(m["pred_eval"]),
                "pred_ignored_border": int(m["pred_ignored_border"]),
                "tp": int(m["tp"]),
                "fp": int(m["fp"]),
                "fn": int(m["fn"]),
                "precision": float(m["precision"]),
                "recall": float(m["recall"]),
                "f1": float(m["f1"]),
                "overseg_ratio": float(m["overseg_ratio_pred_eval_gt"]),
            })

df = pd.DataFrame(rows)

all_csv = OUT / "manual_clean_multi_iou_all_metrics.csv"
df.to_csv(all_csv, index=False)

summary = (
    df.groupby(["method", "threshold"], as_index=False)
    .agg(
        mean_precision=("precision", "mean"),
        mean_recall=("recall", "mean"),
        mean_f1=("f1", "mean"),
        total_tp=("tp", "sum"),
        total_fp=("fp", "sum"),
        total_fn=("fn", "sum"),
        mean_overseg=("overseg_ratio", "mean"),
        n=("f1", "count"),
    )
)

summary_csv = OUT / "manual_clean_multi_iou_summary.csv"
summary.to_csv(summary_csv, index=False)

methods = ["Min-tree-core+WS", "Otsu+WS", "BA-TOS"]
order = {m: i for i, m in enumerate(methods)}

# F1 multi-IoU table
pivot_f1 = summary.pivot(index="method", columns="threshold", values="mean_f1").reset_index()
pivot_f1 = pivot_f1.rename(columns={
    0.25: "F1@0.25",
    0.5: "F1@0.50",
    0.50: "F1@0.50",
    0.75: "F1@0.75",
})
pivot_f1["order"] = pivot_f1["method"].map(order)
pivot_f1 = pivot_f1.sort_values("order").drop(columns=["order"])

pivot_csv = OUT / "table_manual_clean_multi_iou_f1.csv"
pivot_f1.to_csv(pivot_csv, index=False)

# Combined table: compact enough to replace current main table if needed
s050 = summary[np.isclose(summary["threshold"], 0.50)].copy()
s050["order"] = s050["method"].map(order)
s050 = s050.sort_values("order").drop(columns=["order"])

combined = s050[[
    "method",
    "mean_precision",
    "mean_recall",
    "mean_f1",
    "mean_overseg",
    "total_fp",
    "total_fn",
]].copy()

combined = combined.merge(pivot_f1, on="method", how="left")

combined = combined[[
    "method",
    "mean_precision",
    "mean_recall",
    "F1@0.25",
    "F1@0.50",
    "F1@0.75",
    "mean_overseg",
    "total_fp",
    "total_fn",
]]

combined_csv = OUT / "table_manual_clean_combined_results.csv"
combined.to_csv(combined_csv, index=False)

def fmt(x):
    if pd.isna(x):
        return "--"
    return f"{float(x):.3f}"

# LaTeX multi-IoU
latex_multi = []
latex_multi.append(r"\begin{table}[t]")
latex_multi.append(r"\centering")
latex_multi.append(r"\caption{Object-level F1 under multiple IoU thresholds across the twelve evaluated 3D crops.}")
latex_multi.append(r"\label{tab:multi_iou}")
latex_multi.append(r"\begin{tabular}{lrrr}")
latex_multi.append(r"\hline")
latex_multi.append(r"Method & F1@0.25 & F1@0.50 & F1@0.75 \\")
latex_multi.append(r"\hline")
for _, r in pivot_f1.iterrows():
    latex_multi.append(
        f"{r['method']} & {fmt(r.get('F1@0.25'))} & {fmt(r.get('F1@0.50'))} & {fmt(r.get('F1@0.75'))} \\\\"
    )
latex_multi.append(r"\hline")
latex_multi.append(r"\end{tabular}")
latex_multi.append(r"\end{table}")

multi_tex = OUT / "table_manual_clean_multi_iou_f1_latex.tex"
multi_tex.write_text("\n".join(latex_multi), encoding="utf-8")

# LaTeX combined table
latex_combined = []
latex_combined.append(r"\begin{table}[t]")
latex_combined.append(r"\centering")
latex_combined.append(r"\caption{Instance-level comparison across two consecutive scans and twelve representative 3D crops.}")
latex_combined.append(r"\label{tab:main_results_iou050}")
latex_combined.append(r"\resizebox{\columnwidth}{!}{%")
latex_combined.append(r"\begin{tabular}{lrrrrrrrr}")
latex_combined.append(r"\hline")
latex_combined.append(r"Method & P@.50 & R@.50 & F1@.25 & F1@.50 & F1@.75 & Overseg. & FP & FN \\")
latex_combined.append(r"\hline")
for _, r in combined.iterrows():
    latex_combined.append(
        f"{r['method']} & {fmt(r['mean_precision'])} & {fmt(r['mean_recall'])} & {fmt(r['F1@0.25'])} & {fmt(r['F1@0.50'])} & {fmt(r['F1@0.75'])} & {fmt(r['mean_overseg'])} & {int(r['total_fp'])} & {int(r['total_fn'])} \\\\"
    )
latex_combined.append(r"\hline")
latex_combined.append(r"\end{tabular}%")
latex_combined.append(r"}")
latex_combined.append(r"\end{table}")

combined_tex = OUT / "table_manual_clean_combined_results_latex.tex"
combined_tex.write_text("\n".join(latex_combined), encoding="utf-8")

print("\nwrote:", all_csv)
print("wrote:", summary_csv)
print("wrote:", pivot_csv)
print("wrote:", combined_csv)
print("wrote:", multi_tex)
print("wrote:", combined_tex)

print("\n=== MULTI-IOU F1 ===")
print(pivot_f1.to_string(index=False))

print("\n=== COMBINED TABLE ===")
print(combined.to_string(index=False))
