#!/usr/bin/env python3
import csv
from pathlib import Path
import statistics as st

batos_path = Path("results_batos/paper_tables/batos_instance_selected_crops/batos_instance_selected_crops_summary.csv")
otsu_path = Path("results_batos/baselines/otsu_selected_crops/otsu_selected_crops_summary.csv")

out_dir = Path("results_batos/paper_tables")
out_dir.mkdir(parents=True, exist_ok=True)

def read_rows(path, method):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            r["method"] = method
            for k in ["gt_valid", "pred_labels", "fp050", "fn050"]:
                r[k] = int(float(r[k]))
            for k in ["fg_dice", "rec050_bestmatch", "precision050", "recall050", "f1050", "overseg_ratio"]:
                r[k] = float(r[k])
            rows.append(r)
    return rows

rows = read_rows(batos_path, "BA-TOS") + read_rows(otsu_path, "Otsu+WS")

detail_out = out_dir / "batos_vs_otsu_selected_crops.csv"
with open(detail_out, "w", newline="") as f:
    fieldnames = [
        "method", "crop_id", "gt_valid", "pred_labels", "fg_dice",
        "rec050_bestmatch", "precision050", "recall050", "f1050",
        "fp050", "fn050", "overseg_ratio"
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in rows:
        w.writerow({k: r[k] for k in fieldnames})

stats = []
for method in ["BA-TOS", "Otsu+WS"]:
    subset = [r for r in rows if r["method"] == method]
    for metric in ["fg_dice", "precision050", "recall050", "f1050", "overseg_ratio"]:
        vals = [r[metric] for r in subset]
        stats.append({
            "method": method,
            "metric": metric,
            "mean": st.mean(vals),
            "median": st.median(vals),
            "min": min(vals),
            "max": max(vals),
            "std": st.stdev(vals) if len(vals) > 1 else 0.0,
        })

stats_out = out_dir / "batos_vs_otsu_selected_crops_stats.csv"
with open(stats_out, "w", newline="") as f:
    fieldnames = ["method", "metric", "mean", "median", "min", "max", "std"]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in stats:
        w.writerow({
            "method": r["method"],
            "metric": r["metric"],
            "mean": f"{r['mean']:.6f}",
            "median": f"{r['median']:.6f}",
            "min": f"{r['min']:.6f}",
            "max": f"{r['max']:.6f}",
            "std": f"{r['std']:.6f}",
        })

tex_out = out_dir / "table_batos_vs_otsu_stats_latex.tex"
with open(tex_out, "w") as f:
    f.write("\\begin{table}[t]\n")
    f.write("\\centering\n")
    f.write("\\caption{Mean instance-level comparison across six representative EFRGP01\\_00 crops.}\n")
    f.write("\\label{tab:batos_vs_otsu}\n")
    f.write("\\resizebox{\\columnwidth}{!}{%\n")
    f.write("\\begin{tabular}{lrrrrr}\n")
    f.write("\\hline\n")
    f.write("Method & FG Dice & P@0.50 & R@0.50 & F1@0.50 & Overseg. \\\\\n")
    f.write("\\hline\n")

    for method in ["Otsu+WS", "BA-TOS"]:
        d = {r["metric"]: r for r in stats if r["method"] == method}
        f.write(
            f"{method} & "
            f"{d['fg_dice']['mean']:.3f} & "
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
