#!/usr/bin/env python3
from pathlib import Path
import csv
import subprocess
import sys
import pandas as pd
import numpy as np

MANIFEST = Path("results_batos/paper_tables/final_revision/revision_manifest.csv")
OUT = Path("results_batos/paper_tables/final_revision")
METRICS_DIR = OUT / "metrics_runs"

OUT.mkdir(parents=True, exist_ok=True)
METRICS_DIR.mkdir(parents=True, exist_ok=True)

INSTANCE_SCRIPT = Path("code/segmentation/batos/batos_instance_metrics.py")
if not INSTANCE_SCRIPT.exists():
    raise SystemExit(f"missing {INSTANCE_SCRIPT}")

if not MANIFEST.exists():
    raise SystemExit(f"missing {MANIFEST}. Run build_revision_manifest.py first.")

thresholds = ["0.25", "0.50", "0.75"]
rows = []

with MANIFEST.open() as f:
    manifest_rows = list(csv.DictReader(f))

valid_rows = [
    r for r in manifest_rows
    if r["gt_found"] == "1" and r["pred_found"] == "1"
]

print(f"manifest rows: {len(manifest_rows)}")
print(f"valid rows   : {len(valid_rows)}")

for r in valid_rows:
    scan = r["scan"]
    crop_id = r["crop_id"]
    method = r["method"]

    safe_method = (
        method.replace("+", "plus")
              .replace("/", "_")
              .replace(" ", "_")
              .replace("-", "_")
    )

    for th in thresholds:
        safe_th = th.replace(".", "p")
        run_dir = METRICS_DIR / scan / crop_id / safe_method / f"iou_{safe_th}"
        run_dir.mkdir(parents=True, exist_ok=True)

        name = f"{scan}_{crop_id}_{safe_method}_iou{safe_th}"

        cmd = [
            sys.executable, str(INSTANCE_SCRIPT),
            "--gt", r["gt"],
            "--pred", r["pred"],
            "--out-dir", str(run_dir),
            "--min-gt-area", "1000",
            "--ignore-gt-border",
            "--thresholds", th,
            "--name", name,
        ]

        print("RUN:", " ".join(cmd))
        subprocess.run(cmd, check=True)

        metric_file = run_dir / f"{name}_instance_metrics.csv"
        if not metric_file.exists():
            print("missing metric file:", metric_file)
            continue

        df = pd.read_csv(metric_file)

        for _, m in df.iterrows():
            rows.append({
                "scan": scan,
                "crop_id": crop_id,
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

if not rows:
    raise SystemExit("No metric rows generated. Check manifest paths.")

all_df = pd.DataFrame(rows)
all_csv = OUT / "revision_all_instance_metrics.csv"
all_df.to_csv(all_csv, index=False)
print("wrote:", all_csv)

summary = (
    all_df
    .groupby(["method", "threshold"], as_index=False)
    .agg(
        mean_precision=("precision", "mean"),
        mean_recall=("recall", "mean"),
        mean_f1=("f1", "mean"),
        median_f1=("f1", "median"),
        std_f1=("f1", "std"),
        total_tp=("tp", "sum"),
        total_fp=("fp", "sum"),
        total_fn=("fn", "sum"),
        mean_overseg=("overseg_ratio", "mean"),
        n=("f1", "count"),
    )
)

summary_csv = OUT / "revision_summary_by_method_threshold.csv"
summary.to_csv(summary_csv, index=False)
print("wrote:", summary_csv)

main_methods = ["Min-tree-core+WS", "Otsu+WS", "BA-TOS"]

multi = summary[summary["method"].isin(main_methods)].copy()
pivot = multi.pivot(index="method", columns="threshold", values="mean_f1").reset_index()
pivot = pivot.rename(columns={
    0.25: "F1@0.25",
    0.5: "F1@0.50",
    0.50: "F1@0.50",
    0.75: "F1@0.75",
})
order = {m: i for i, m in enumerate(main_methods)}
pivot["order"] = pivot["method"].map(order)
pivot = pivot.sort_values("order").drop(columns=["order"])

multi_csv = OUT / "table_multi_iou_f1.csv"
pivot.to_csv(multi_csv, index=False)
print("wrote:", multi_csv)

primary = summary[
    (summary["method"].isin(main_methods)) &
    (np.isclose(summary["threshold"], 0.50))
].copy()

primary["order"] = primary["method"].map(order)
primary = primary.sort_values("order")

primary_out = primary[[
    "method", "mean_precision", "mean_recall", "mean_f1",
    "mean_overseg", "total_fp", "total_fn", "n"
]]

primary_csv = OUT / "table_primary_results_iou050.csv"
primary_out.to_csv(primary_csv, index=False)
print("wrote:", primary_csv)

abl_methods = ["BA-TOS", "BA-TOS raw/no-size-filter"]
abl = summary[
    (summary["method"].isin(abl_methods)) &
    (np.isclose(summary["threshold"], 0.50))
].copy()

abl_order = {m: i for i, m in enumerate(abl_methods)}
abl["order"] = abl["method"].map(abl_order)
abl = abl.sort_values("order")

abl_out = abl[[
    "method", "mean_precision", "mean_recall", "mean_f1",
    "mean_overseg", "total_fp", "total_fn", "n"
]]

abl_csv = OUT / "table_ablation_size_filter_iou050.csv"
abl_out.to_csv(abl_csv, index=False)
print("wrote:", abl_csv)

def fmt(x):
    if pd.isna(x):
        return "--"
    return f"{x:.3f}"

# LaTeX primary table
latex_primary = []
latex_primary.append(r"\begin{table}[t]")
latex_primary.append(r"\centering")
latex_primary.append(r"\caption{Instance-level comparison across the evaluated 3D crops at IoU@0.50.}")
latex_primary.append(r"\label{tab:main_results_iou050}")
latex_primary.append(r"\resizebox{\columnwidth}{!}{%")
latex_primary.append(r"\begin{tabular}{lrrrrrr}")
latex_primary.append(r"\hline")
latex_primary.append(r"Method & P@0.50 & R@0.50 & F1@0.50 & Overseg. & FP & FN \\")
latex_primary.append(r"\hline")
for _, r in primary_out.iterrows():
    latex_primary.append(
        f"{r['method']} & {fmt(r['mean_precision'])} & {fmt(r['mean_recall'])} & {fmt(r['mean_f1'])} & {fmt(r['mean_overseg'])} & {int(r['total_fp'])} & {int(r['total_fn'])} \\\\"
    )
latex_primary.append(r"\hline")
latex_primary.append(r"\end{tabular}%")
latex_primary.append(r"}")
latex_primary.append(r"\end{table}")

(OUT / "table_primary_results_iou050_latex.tex").write_text(
    "\n".join(latex_primary), encoding="utf-8"
)
print("wrote:", OUT / "table_primary_results_iou050_latex.tex")

# LaTeX multi IoU
latex_multi = []
latex_multi.append(r"\begin{table}[t]")
latex_multi.append(r"\centering")
latex_multi.append(r"\caption{Object-level F1 under multiple IoU thresholds.}")
latex_multi.append(r"\label{tab:multi_iou}")
latex_multi.append(r"\begin{tabular}{lrrr}")
latex_multi.append(r"\hline")
latex_multi.append(r"Method & F1@0.25 & F1@0.50 & F1@0.75 \\")
latex_multi.append(r"\hline")
for _, r in pivot.iterrows():
    latex_multi.append(
        f"{r['method']} & {fmt(r.get('F1@0.25'))} & {fmt(r.get('F1@0.50'))} & {fmt(r.get('F1@0.75'))} \\\\"
    )
latex_multi.append(r"\hline")
latex_multi.append(r"\end{tabular}")
latex_multi.append(r"\end{table}")

(OUT / "table_multi_iou_f1_latex.tex").write_text(
    "\n".join(latex_multi), encoding="utf-8"
)
print("wrote:", OUT / "table_multi_iou_f1_latex.tex")

# LaTeX ablation
latex_abl = []
latex_abl.append(r"\begin{table}[t]")
latex_abl.append(r"\centering")
latex_abl.append(r"\caption{Initial ablation of the size-aware post-filter at IoU@0.50.}")
latex_abl.append(r"\label{tab:ablation_size_filter}")
latex_abl.append(r"\resizebox{\columnwidth}{!}{%")
latex_abl.append(r"\begin{tabular}{lrrrrrr}")
latex_abl.append(r"\hline")
latex_abl.append(r"Variant & P@0.50 & R@0.50 & F1@0.50 & Overseg. & FP & FN \\")
latex_abl.append(r"\hline")
for _, r in abl_out.iterrows():
    name = "Full BA-TOS" if r["method"] == "BA-TOS" else "w/o size filter"
    latex_abl.append(
        f"{name} & {fmt(r['mean_precision'])} & {fmt(r['mean_recall'])} & {fmt(r['mean_f1'])} & {fmt(r['mean_overseg'])} & {int(r['total_fp'])} & {int(r['total_fn'])} \\\\"
    )
latex_abl.append(r"\hline")
latex_abl.append(r"\end{tabular}%")
latex_abl.append(r"}")
latex_abl.append(r"\end{table}")

(OUT / "table_ablation_size_filter_iou050_latex.tex").write_text(
    "\n".join(latex_abl), encoding="utf-8"
)
print("wrote:", OUT / "table_ablation_size_filter_iou050_latex.tex")

print()
print("=== PRIMARY @0.50 ===")
print(primary_out.to_string(index=False))

print()
print("=== MULTI IOU ===")
print(pivot.to_string(index=False))

print()
print("=== ABLATION SIZE FILTER ===")
print(abl_out.to_string(index=False))
