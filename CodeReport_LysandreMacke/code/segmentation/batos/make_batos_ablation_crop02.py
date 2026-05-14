#!/usr/bin/env python3
import csv
from pathlib import Path
import sys


ITEMS = [
    (
        "BA-TOS markers + Macke-style watershed",
        Path("results_batos/gustavo_EFR_00_crop_x600_y800_z300/eval_interior_min1000/summary.csv"),
    ),
    (
        "BA-TOS markers + distance watershed",
        Path("results_batos/gustavo_EFR_00_crop_x600_y800_z300/eval_batos_pipeline_final/summary.csv"),
    ),
    (
        "BA-TOS markers + GT-mask watershed oracle",
        Path("results_batos/gustavo_EFR_00_crop_x600_y800_z300/eval_oracle_gtmask_interior_min1000/summary.csv"),
    ),
]


def read_summary(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing summary file: {path}")

    vals = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            vals[row["metric"]] = row["value"]
    return vals


def fnum(vals, key, default="0"):
    try:
        return float(vals.get(key, default))
    except Exception:
        return 0.0


def main():
    out_dir = Path("results_batos/paper_tables")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []

    for method, path in ITEMS:
        vals = read_summary(path)

        rows.append({
            "method": method,
            "gt_valid": int(fnum(vals, "gt_valid_ids_after_filter")),
            "pred_labels": int(fnum(vals, "pred_label_count")),
            "fg_iou": fnum(vals, "foreground_iou"),
            "fg_dice": fnum(vals, "foreground_dice"),
            "fg_precision": fnum(vals, "foreground_precision"),
            "fg_recall": fnum(vals, "foreground_recall"),
            "rec025": fnum(vals, "gt_object_recall_iou_0.25"),
            "rec050": fnum(vals, "gt_object_recall_iou_0.5"),
            "matched050": int(fnum(vals, "gt_matched_count_iou_0.5")),
            "source": str(path),
        })

    csv_out = out_dir / "batos_ablation_crop02.csv"

    with open(csv_out, "w", newline="") as f:
        fieldnames = [
            "method",
            "gt_valid",
            "pred_labels",
            "fg_iou",
            "fg_dice",
            "fg_precision",
            "fg_recall",
            "rec025",
            "rec050",
            "matched050",
            "source",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    tex_out = out_dir / "table_batos_ablation_crop02_latex.tex"

    with open(tex_out, "w") as f:
        f.write("\\begin{table}[t]\n")
        f.write("\\centering\n")
        f.write("\\caption{Ablation analysis on crop02. The GT-mask oracle is reported only as a diagnostic upper bound.}\n")
        f.write("\\label{tab:batos_ablation_crop02}\n")
        f.write("\\resizebox{\\columnwidth}{!}{%\n")
        f.write("\\begin{tabular}{lrrrrr}\n")
        f.write("\\hline\n")
        f.write("Method & GT & Pred. & FG Dice & Rec.@0.25 & Rec.@0.50 \\\\\n")
        f.write("\\hline\n")

        for r in rows:
            f.write(
                f"{r['method']} & "
                f"{r['gt_valid']} & "
                f"{r['pred_labels']} & "
                f"{r['fg_dice']:.3f} & "
                f"{r['rec025']:.3f} & "
                f"{r['rec050']:.3f} \\\\\n"
            )

        f.write("\\hline\n")
        f.write("\\end{tabular}%\n")
        f.write("}\n")
        f.write("\\end{table}\n")

    print(f"wrote: {csv_out}")
    print(f"wrote: {tex_out}")

    print("")
    print("=== Ablation summary ===")
    for r in rows:
        print(
            f"{r['method']}: "
            f"GT={r['gt_valid']}, "
            f"Pred={r['pred_labels']}, "
            f"Dice={r['fg_dice']:.4f}, "
            f"Rec@0.25={r['rec025']:.4f}, "
            f"Rec@0.50={r['rec050']:.4f}, "
            f"Matched@0.50={r['matched050']}"
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("[ERROR]", e, file=sys.stderr)
        sys.exit(1)
