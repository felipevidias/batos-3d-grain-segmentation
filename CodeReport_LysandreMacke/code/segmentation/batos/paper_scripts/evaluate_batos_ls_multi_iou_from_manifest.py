#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys
import pandas as pd

ROOT = Path(".")
MANIFEST = ROOT / "results_batos/paper_tables/batos_local_split_large_labels/local_split_manifest.csv"
METRICS = ROOT / "code/segmentation/batos/batos_instance_metrics.py"

OUT = ROOT / "results_batos/paper_tables/batos_local_split_large_labels"
OUT.mkdir(parents=True, exist_ok=True)

THRESHOLDS = ["0.25", "0.50", "0.75"]


def require(path):
    if not path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")


def get_value(row, candidates):
    for c in candidates:
        if c in row:
            return row[c]
    raise KeyError(f"Nenhuma coluna encontrada entre: {candidates}. Colunas disponíveis: {list(row.keys())}")


def run_metric(gt, pred, out_dir, name, threshold):
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, str(METRICS),
        "--gt", str(gt),
        "--pred", str(pred),
        "--out-dir", str(out_dir),
        "--min-gt-area", "1000",
        "--ignore-gt-border",
        "--thresholds", threshold,
        "--name", name,
    ]

    print("RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    metric_file = out_dir / f"{name}_instance_metrics.csv"

    if not metric_file.exists():
        raise FileNotFoundError(f"Métrica não gerada: {metric_file}")

    df = pd.read_csv(metric_file)

    if len(df) == 0:
        raise RuntimeError(f"CSV de métrica vazio: {metric_file}")

    return df.iloc[0].to_dict()


def main():
    require(MANIFEST)
    require(METRICS)

    manifest = pd.read_csv(MANIFEST)

    all_rows = []

    for _, r in manifest.iterrows():
        scan = str(r["scan"])
        crop_id = str(r["crop_id"])
        gt = Path(str(r["gt"]))
        pred = Path(str(r["split_pred"]))

        require(gt)
        require(pred)

        for th in THRESHOLDS:
            tag = th.replace(".", "p")
            name = f"{scan}_{crop_id}_BA_TOS_LS_iou_{tag}"
            metric_out = OUT / scan / crop_id / f"metrics_iou_{tag}"

            m = run_metric(
                gt=gt,
                pred=pred,
                out_dir=metric_out,
                name=name,
                threshold=th,
            )

            row = {
                "method": "BA-TOS-LS",
                "scan": scan,
                "crop_id": crop_id,
                "threshold": float(th),
                "gt": str(gt),
                "pred": str(pred),
                **m,
            }

            all_rows.append(row)

    metrics_df = pd.DataFrame(all_rows)
    metrics_csv = OUT / "local_split_multi_iou_crop_metrics.csv"
    metrics_df.to_csv(metrics_csv, index=False)

    summary_rows = []

    for th in [0.25, 0.50, 0.75]:
        sub = metrics_df[metrics_df["threshold"] == th].copy()

        if len(sub) == 0:
            raise RuntimeError(f"Nenhuma métrica encontrada para threshold={th}")

        summary_rows.append({
            "method": "BA-TOS-LS",
            "threshold": th,
            "mean_precision": sub["precision"].astype(float).mean(),
            "mean_recall": sub["recall"].astype(float).mean(),
            "mean_f1": sub["f1"].astype(float).mean(),
            "mean_overseg": sub["overseg_ratio_pred_eval_gt"].astype(float).mean(),
            "total_fp": int(sub["fp"].astype(int).sum()),
            "total_fn": int(sub["fn"].astype(int).sum()),
            "n_crops": int(len(sub)),
        })

    by_threshold = pd.DataFrame(summary_rows)
    by_threshold_csv = OUT / "local_split_multi_iou_by_threshold.csv"
    by_threshold.to_csv(by_threshold_csv, index=False)

    row025 = by_threshold[by_threshold["threshold"] == 0.25].iloc[0]
    row050 = by_threshold[by_threshold["threshold"] == 0.50].iloc[0]
    row075 = by_threshold[by_threshold["threshold"] == 0.75].iloc[0]

    final = pd.DataFrame([{
        "method": "BA-TOS-LS",
        "mean_precision": row050["mean_precision"],
        "mean_recall": row050["mean_recall"],
        "F1@0.25": row025["mean_f1"],
        "F1@0.50": row050["mean_f1"],
        "F1@0.75": row075["mean_f1"],
        "mean_overseg": row050["mean_overseg"],
        "total_fp": int(row050["total_fp"]),
        "total_fn": int(row050["total_fn"]),
        "n_crops": int(row050["n_crops"]),
    }])

    final_csv = OUT / "local_split_multi_iou_summary.csv"
    final.to_csv(final_csv, index=False)

    print()
    print("=== BA-TOS-LS multi-IoU summary ===")
    print(final.to_string(index=False))
    print()
    print("wrote:", metrics_csv)
    print("wrote:", by_threshold_csv)
    print("wrote:", final_csv)


if __name__ == "__main__":
    main()
