#!/usr/bin/env python3
import csv
from pathlib import Path
import statistics as st

batos_sweep = Path("results_batos/paper_tables/batos_size_filter_sweep/batos_size_filter_sweep_summary.csv")
otsu_path = Path("results_batos/baselines/otsu_selected_crops/otsu_selected_crops_summary.csv")

out_dir = Path("results_batos/paper_tables")
out_dir.mkdir(parents=True, exist_ok=True)

rows = []

with open(batos_sweep) as f:
    for r in csv.DictReader(f):
        if int(r["min_size"]) != 500:
            continue
        rows.append({
            "method": "BA-TOS",
            "crop_id": r["crop_id"],
            "gt_valid": int(r["gt_valid"]),
            "pred_labels": int(r["pred_total"]),
            "fg_dice": "",
            "precision050": float(r["precision050"]),
            "recall050": float(r["recall050"]),
            "f1050": float(r["f1050"]),
            "fp050": int(r["fp050"]),
            "fn050": int(r["fn050"]),
            "overseg_ratio": float(r["overseg_ratio"]),
        })

with open(otsu_path) as f:
    for r in csv.DictReader(f):
        rows.append({
            "method": "Otsu+WS",
            "crop_id": r["crop_id"],
            "gt_valid": int(r["gt_valid"]),
            "pred_labels": int(r["pred_labels"]),
            "fg_dice": float(r["fg_dice"]),
            "precision050": float(r["precision050"]),
            "recall050": float(r["recall050"]),
            "f1050": float(r["f1050"]),
            "fp050": int(r["fp050"]),
            "fn050": int(r["fn050"]),
            "overseg_ratio": float(r["overseg_ratio"]),
        })

detail_out = out_dir / "batos_min500_vs_otsu_selected_crops.csv"
with open(detail_out, "w", newline="") as f:
    fieldnames = [
        "method", "crop_id", "gt_valid", "pred_labels",
        "precision050", "recall050", "f1050",
        "fp050", "fn050", "overseg_ratio"
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in rows:
        w.writerow({k: r[k] for k in fieldnames})

stats_rows = []
for method in ["Otsu+WS", "BA-TOS"]:
    subset = [r for r in rows if r["method"] == method]
    for metric in ["precision050", "recall050", "f1050", "overseg_ratio"]:
        vals = [r[metric] for r in subset]
        stats_rows.append({
            "method": method,
            "metric": metric,
            "mean": st.mean(vals),
            "median": st.median(vals),
            "min": min(vals),
            "max": max(vals),
            "std": st.stdev(vals) if len(vals) > 1 else 0.0,
        })

stats_out = out_dir / "batos_min500_vs_otsu_stats.csv"
with open(stats_out, "w", newline="") as f:
    fieldnames = ["method", "metric", "mean", "median", "min", "max", "std"]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in stats_rows:
        w.writerow({
            "method": r["method"],
            "metric": r["metric"],
            "mean": f"{r['mean']:.6f}",
            "median": f"{r['median']:.6f}",
            "min": f"{r['min']:.6f}",
            "max": f"{r['max']:.6f}",
            "std": f"{r['std']:.6f}",
        })

tex_out = out_dir / "table_batos_min500_vs_otsu_latex.tex"
with open(tex_out, "w") as f:
    f.write("\\begin{table}[t]\n")
    f.write("\\centering\n")
    f.write("\\caption{Instance-level comparison across six representative EFRGP01\\_00 crops. BA-TOS includes size-aware post-filtering with a minimum predicted-instance size of 500 voxels.}\n")
    f.write("\\label{tab:batos_min500_vs_otsu}\n")
    f.write("\\resizebox{\\columnwidth}{!}{%\n")
    f.write("\\begin{tabular}{lrrrr}\n")
    f.write("\\hline\n")
    f.write("Method & P@0.50 & R@0.50 & F1@0.50 & Overseg. \\\\\n")
    f.write("\\hline\n")

    for method in ["Otsu+WS", "BA-TOS"]:
        d = {r["metric"]: r for r in stats_rows if r["method"] == method}
        f.write(
            f"{method} & "
            f"{d['precision050']['mean']:.3f} & "
            f"{d['recall050']['mean']:.3f} & "
            f"{d['f1050']['mean']:.3f} & "
            f"{d['overseg_ratio']['mean']:.3f} \\\\\n"
        )

    f.write("\\hline\n")
    f.write("\\end{tabular}%\n")
    f.write("}\n")
    f.write("\\end{table}\n")

print("wrote:", detail_out)
print("wrote:", stats_out)
print("wrote:", tex_out)
