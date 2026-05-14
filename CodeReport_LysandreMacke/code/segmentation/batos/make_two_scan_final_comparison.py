#!/usr/bin/env python3
import csv
from pathlib import Path
import statistics as st

out_dir = Path("results_batos/paper_tables")
out_dir.mkdir(parents=True, exist_ok=True)

rows = []

# Scan 1: EFRGP01_00
p00 = Path("results_batos/paper_tables/final_baseline_comparison_selected_crops.csv")
with open(p00) as f:
    for r in csv.DictReader(f):
        rows.append({
            "scan": "EFRGP01_00",
            "crop_id": r["crop_id"],
            "method": r["method"],
            "gt_valid": int(float(r["gt_valid"])),
            "pred_labels": int(float(r["pred_labels"])),
            "precision050": float(r["precision050"]),
            "recall050": float(r["recall050"]),
            "f1050": float(r["f1050"]),
            "fp050": int(float(r["fp050"])),
            "fn050": int(float(r["fn050"])),
            "overseg_ratio": float(r["overseg_ratio"]),
        })

# Scan 2: EFRGP01_01
p01 = Path("results_batos/EFRGP01_01_selected_crops/EFRGP01_01_selected_crops_methods_summary.csv")
with open(p01) as f:
    for r in csv.DictReader(f):
        rows.append({
            "scan": r["scan"],
            "crop_id": r["crop_id"],
            "method": r["method"],
            "gt_valid": int(float(r["gt_valid"])),
            "pred_labels": int(float(r["pred_labels"])),
            "precision050": float(r["precision050"]),
            "recall050": float(r["recall050"]),
            "f1050": float(r["f1050"]),
            "fp050": int(float(r["fp050"])),
            "fn050": int(float(r["fn050"])),
            "overseg_ratio": float(r["overseg_ratio"]),
        })

detail_out = out_dir / "two_scan_final_comparison_all_crops.csv"
with open(detail_out, "w", newline="") as f:
    fieldnames = [
        "scan", "crop_id", "method", "gt_valid", "pred_labels",
        "precision050", "recall050", "f1050",
        "fp050", "fn050", "overseg_ratio"
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in rows:
        w.writerow({k: r[k] for k in fieldnames})

stats_rows = []

for group_name, group_rows in [
    ("ALL", rows),
    ("EFRGP01_00", [r for r in rows if r["scan"] == "EFRGP01_00"]),
    ("EFRGP01_01", [r for r in rows if r["scan"] == "EFRGP01_01"]),
]:
    for method in ["Macke-MinTree+WS", "Otsu+WS", "BA-TOS"]:
        subset = [r for r in group_rows if r["method"] == method]
        if not subset:
            continue

        for metric in ["precision050", "recall050", "f1050", "overseg_ratio"]:
            vals = [r[metric] for r in subset]
            stats_rows.append({
                "group": group_name,
                "method": method,
                "metric": metric,
                "mean": st.mean(vals),
                "median": st.median(vals),
                "min": min(vals),
                "max": max(vals),
                "std": st.stdev(vals) if len(vals) > 1 else 0.0,
            })

        stats_rows.append({
            "group": group_name,
            "method": method,
            "metric": "total_fp",
            "mean": sum(r["fp050"] for r in subset),
            "median": "",
            "min": "",
            "max": "",
            "std": "",
        })

        stats_rows.append({
            "group": group_name,
            "method": method,
            "metric": "total_fn",
            "mean": sum(r["fn050"] for r in subset),
            "median": "",
            "min": "",
            "max": "",
            "std": "",
        })

stats_out = out_dir / "two_scan_final_comparison_stats.csv"
with open(stats_out, "w", newline="") as f:
    fieldnames = ["group", "method", "metric", "mean", "median", "min", "max", "std"]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in stats_rows:
        out = dict(r)
        for k in ["mean", "median", "min", "max", "std"]:
            if isinstance(out[k], float):
                out[k] = f"{out[k]:.6f}"
        w.writerow(out)

tex_out = out_dir / "table_two_scan_final_comparison_latex.tex"
with open(tex_out, "w") as f:
    f.write("\\begin{table}[t]\n")
    f.write("\\centering\n")
    f.write("\\caption{Instance-level comparison across two consecutive scans and twelve representative 3D crops.}\n")
    f.write("\\label{tab:two_scan_final_comparison}\n")
    f.write("\\resizebox{\\columnwidth}{!}{%\n")
    f.write("\\begin{tabular}{lrrrrrr}\n")
    f.write("\\hline\n")
    f.write("Method & P@0.50 & R@0.50 & F1@0.50 & Overseg. & FP & FN \\\\\n")
    f.write("\\hline\n")

    for method in ["Macke-MinTree+WS", "Otsu+WS", "BA-TOS"]:
        d = {
            r["metric"]: r
            for r in stats_rows
            if r["group"] == "ALL" and r["method"] == method
        }
        f.write(
            f"{method} & "
            f"{float(d['precision050']['mean']):.3f} & "
            f"{float(d['recall050']['mean']):.3f} & "
            f"{float(d['f1050']['mean']):.3f} & "
            f"{float(d['overseg_ratio']['mean']):.3f} & "
            f"{int(d['total_fp']['mean'])} & "
            f"{int(d['total_fn']['mean'])} \\\\\n"
        )

    f.write("\\hline\n")
    f.write("\\end{tabular}%\n")
    f.write("}\n")
    f.write("\\end{table}\n")

print("wrote:", detail_out)
print("wrote:", stats_out)
print("wrote:", tex_out)

print("\n=== Overall summary ===")
for method in ["Macke-MinTree+WS", "Otsu+WS", "BA-TOS"]:
    d = {
        r["metric"]: r
        for r in stats_rows
        if r["group"] == "ALL" and r["method"] == method
    }
    print(
        method,
        "P=", f"{float(d['precision050']['mean']):.4f}",
        "R=", f"{float(d['recall050']['mean']):.4f}",
        "F1=", f"{float(d['f1050']['mean']):.4f}",
        "Over=", f"{float(d['overseg_ratio']['mean']):.4f}",
        "FP=", int(d["total_fp"]["mean"]),
        "FN=", int(d["total_fn"]["mean"]),
    )
