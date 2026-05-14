#!/usr/bin/env python3
import csv
from pathlib import Path


OUT = Path("results_batos/paper_tables")
OUT.mkdir(parents=True, exist_ok=True)

methods = [
    {
        "method": "Otsu + distance watershed",
        "fg_summary": Path("results_batos/baselines/crop02_otsu_ws_eval/summary.csv"),
        "inst": Path("results_batos/paper_tables/instance_metrics_crop02/otsu_ws_crop02_instance_metrics.csv"),
    },
    {
        "method": "BA-TOS + Macke-style watershed",
        "fg_summary": Path("results_batos/gustavo_EFR_00_crop_x600_y800_z300/eval_interior_min1000/summary.csv"),
        "inst": None,
    },
    {
        "method": "BA-TOS + distance watershed",
        "fg_summary": Path("results_batos/gustavo_EFR_00_crop_x600_y800_z300/eval_batos_pipeline_final/summary.csv"),
        "inst": Path("results_batos/paper_tables/instance_metrics_crop02/batos_distance_crop02_instance_metrics.csv"),
    },
    {
        "method": "GT-mask oracle",
        "fg_summary": Path("results_batos/gustavo_EFR_00_crop_x600_y800_z300/eval_oracle_gtmask_interior_min1000/summary.csv"),
        "inst": None,
    },
]


def read_summary(path):
    vals = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            vals[r["metric"]] = r["value"]
    return vals


def read_inst_050(path):
    if path is None or not path.exists():
        return None

    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if abs(float(r["threshold"]) - 0.5) < 1e-9:
                return r
    return None


rows = []

for m in methods:
    s = read_summary(m["fg_summary"])
    inst = read_inst_050(m["inst"])

    row = {
        "method": m["method"],
        "gt_valid": int(float(s.get("gt_valid_ids_after_filter", 0))),
        "pred_labels": int(float(s.get("pred_label_count", 0))),
        "fg_dice": float(s.get("foreground_dice", 0)),
        "rec050_bestmatch": float(s.get("gt_object_recall_iou_0.5", 0)),
        "matched050_bestmatch": int(float(s.get("gt_matched_count_iou_0.5", 0))),
        "precision050_hungarian": "",
        "recall050_hungarian": "",
        "f1050_hungarian": "",
        "fp050": "",
        "fn050": "",
        "overseg_ratio": "",
    }

    if inst is not None:
        row.update({
            "precision050_hungarian": float(inst["precision"]),
            "recall050_hungarian": float(inst["recall"]),
            "f1050_hungarian": float(inst["f1"]),
            "fp050": int(inst["fp"]),
            "fn050": int(inst["fn"]),
            "overseg_ratio": float(inst["overseg_ratio_pred_eval_gt"]),
        })

    rows.append(row)


csv_out = OUT / "crop02_method_comparison.csv"
with open(csv_out, "w", newline="") as f:
    fieldnames = [
        "method",
        "gt_valid",
        "pred_labels",
        "fg_dice",
        "rec050_bestmatch",
        "matched050_bestmatch",
        "precision050_hungarian",
        "recall050_hungarian",
        "f1050_hungarian",
        "fp050",
        "fn050",
        "overseg_ratio",
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)


tex_out = OUT / "table_crop02_method_comparison_latex.tex"
with open(tex_out, "w") as f:
    f.write("\\begin{table}[t]\n")
    f.write("\\centering\n")
    f.write("\\caption{Crop02 comparison against reference particle labels. The GT-mask oracle is shown only as a diagnostic upper bound.}\n")
    f.write("\\label{tab:crop02_comparison}\n")
    f.write("\\resizebox{\\columnwidth}{!}{%\n")
    f.write("\\begin{tabular}{lrrrrrr}\n")
    f.write("\\hline\n")
    f.write("Method & Pred. & FG Dice & Rec.@0.50 & P@0.50 & F1@0.50 & Overseg. \\\\\n")
    f.write("\\hline\n")

    for r in rows:
        p = "--" if r["precision050_hungarian"] == "" else f"{r['precision050_hungarian']:.3f}"
        f1 = "--" if r["f1050_hungarian"] == "" else f"{r['f1050_hungarian']:.3f}"
        over = "--" if r["overseg_ratio"] == "" else f"{r['overseg_ratio']:.3f}"

        f.write(
            f"{r['method']} & "
            f"{r['pred_labels']} & "
            f"{r['fg_dice']:.3f} & "
            f"{r['rec050_bestmatch']:.3f} & "
            f"{p} & "
            f"{f1} & "
            f"{over} \\\\\n"
        )

    f.write("\\hline\n")
    f.write("\\end{tabular}%\n")
    f.write("}\n")
    f.write("\\end{table}\n")


print("wrote:", csv_out)
print("wrote:", tex_out)
