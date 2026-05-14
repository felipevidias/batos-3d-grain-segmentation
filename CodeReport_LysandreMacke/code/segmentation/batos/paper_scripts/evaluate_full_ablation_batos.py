#!/usr/bin/env python3
from pathlib import Path
import csv
import subprocess
import sys
import pandas as pd

ROOT = Path(".")
TABLEDIR = ROOT / "results_batos/paper_tables/final_revision_clean"
MANIFEST = TABLEDIR / "full_ablation_manifest.csv"
RUNS = TABLEDIR / "full_ablation_metrics_runs"
INSTANCE = ROOT / "code/segmentation/batos/batos_instance_metrics.py"

RUNS.mkdir(parents=True, exist_ok=True)

if not MANIFEST.exists():
    raise SystemExit(f"missing manifest: {MANIFEST}")

if not INSTANCE.exists():
    raise SystemExit(f"missing instance metric script: {INSTANCE}")

rows = []

manifest = list(csv.DictReader(open(MANIFEST)))
valid = [r for r in manifest if r["gt_found"] == "1" and r["pred_found"] == "1"]

print("manifest rows:", len(manifest))
print("valid rows   :", len(valid))

if len(valid) != len(manifest):
    print("Missing entries:")
    for r in manifest:
        if r["gt_found"] != "1" or r["pred_found"] != "1":
            print(r)
    raise SystemExit("manifest has missing files")

for r in valid:
    scan = r["scan"]
    crop = r["crop_id"]
    variant = r["variant"]

    safe_variant = (
        variant.replace(" ", "_")
        .replace("/", "_")
        .replace("+", "plus")
        .replace("-", "_")
    )

    outdir = RUNS / scan / crop / safe_variant
    outdir.mkdir(parents=True, exist_ok=True)

    name = f"{scan}_{crop}_{safe_variant}"

    cmd = [
        sys.executable, str(INSTANCE),
        "--gt", r["gt"],
        "--pred", r["pred"],
        "--out-dir", str(outdir),
        "--min-gt-area", "1000",
        "--ignore-gt-border",
        "--thresholds", "0.50",
        "--name", name,
    ]

    print("RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    mf = outdir / f"{name}_instance_metrics.csv"
    if not mf.exists():
        raise RuntimeError(f"metric file not found: {mf}")

    df = pd.read_csv(mf)
    m = df.iloc[0]

    rows.append({
        "scan": scan,
        "crop_id": crop,
        "variant": variant,
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
        "overseg": float(m["overseg_ratio_pred_eval_gt"]),
    })

df = pd.DataFrame(rows)

all_csv = TABLEDIR / "full_ablation_all_crops.csv"
df.to_csv(all_csv, index=False)

summary = (
    df.groupby("variant", as_index=False)
    .agg(
        mean_precision=("precision", "mean"),
        mean_recall=("recall", "mean"),
        mean_f1=("f1", "mean"),
        std_f1=("f1", "std"),
        mean_overseg=("overseg", "mean"),
        total_tp=("tp", "sum"),
        total_fp=("fp", "sum"),
        total_fn=("fn", "sum"),
        n=("f1", "count"),
    )
)

order = {
    "Full BA-TOS": 0,
    "w/o boundary rejection": 1,
    "w/o size filter": 2,
    "w/o flattening": 3,
}

summary["order"] = summary["variant"].map(order)
summary = summary.sort_values("order").drop(columns=["order"])

summary_csv = TABLEDIR / "table_full_ablation_batos.csv"
summary.to_csv(summary_csv, index=False)

def fmt(x):
    return f"{float(x):.3f}"

latex = []
latex.append(r"\begin{table}[t]")
latex.append(r"\centering")
latex.append(r"\caption{Ablation study of BA-TOS components at IoU@0.50.}")
latex.append(r"\label{tab:ablation}")
latex.append(r"\resizebox{\columnwidth}{!}{%")
latex.append(r"\begin{tabular}{lrrrrrr}")
latex.append(r"\hline")
latex.append(r"Variant & P@.50 & R@.50 & F1@.50 & Overseg. & FP & FN \\")
latex.append(r"\hline")

for _, r in summary.iterrows():
    latex.append(
        f"{r['variant']} & {fmt(r['mean_precision'])} & {fmt(r['mean_recall'])} & {fmt(r['mean_f1'])} & {fmt(r['mean_overseg'])} & {int(r['total_fp'])} & {int(r['total_fn'])} \\\\"
    )

latex.append(r"\hline")
latex.append(r"\end{tabular}%")
latex.append(r"}")
latex.append(r"\end{table}")

latex_path = TABLEDIR / "table_full_ablation_batos_latex.tex"
latex_path.write_text("\n".join(latex), encoding="utf-8")

print()
print("wrote:", all_csv)
print("wrote:", summary_csv)
print("wrote:", latex_path)

print()
print("=== FULL ABLATION SUMMARY ===")
print(summary.to_string(index=False))
