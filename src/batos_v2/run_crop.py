"""Stage 6: run stages 1-5 on one crop of crops_sibgrapi.json and record provenance.

Writes C:\\BA_TOS_DATA\\runs\\batos_paper_v2\\<crop_id>\\ (01_preprocess, 02_tos, 03_gate, 04_watershed,
05_eval, run.json) and appends one row per (method, size filter, protocol) to results/all_runs.csv.

Methods evaluated:
  batos_tos  - BA-TOS_ToS route (ToS MIN markers)
  batos_sel  - BA-TOS_sel (gate-routed markers; identical to batos_tos when the gate is false)
  edt_otsu   - matched-count EDT-Otsu markers under the same mask (Table II "Otsu+WS" control)
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import evaluate, gate, preprocess, tos_select, watershed
from .common import REPO_ROOT, config_hash, load_config, read3d, write_json

METHOD_MARKERS = {"batos_tos": "tos_e11_min", "batos_sel": "batos_selective", "edt_otsu": "edt_otsu"}


def git_state() -> dict[str, Any]:
    def run(*args):
        return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True).stdout.strip()
    return {"commit": run("rev-parse", "HEAD"), "branch": run("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty": bool(run("status", "--porcelain", "--untracked-files=no"))}


def peak_memory_mb() -> float | None:
    try:
        import psutil
        info = psutil.Process().memory_info()
        return float(getattr(info, "peak_wset", getattr(info, "rss", 0))) / 2 ** 20
    except Exception:
        return None


def find_crop(config: dict[str, Any], crop_id: str) -> dict[str, Any]:
    meta_path = REPO_ROOT / config["paths"]["crops_metadata"]
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if meta.get("crop_convention", {}).get("slicing_status") != "verified":
        raise RuntimeError("crops_sibgrapi.json axis order is not verified")
    for crop in meta["crops"]:
        if crop["crop_id"] == crop_id:
            return crop
    raise KeyError(f"crop_id {crop_id!r} not in {meta_path}")


def append_all_runs(rows: list[dict[str, Any]], config: dict[str, Any]) -> Path:
    path = REPO_ROOT / config["paths"].get("all_runs", "results/all_runs.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    if path.exists():
        frame = pd.concat([pd.read_csv(path), frame], ignore_index=True)
    frame.to_csv(path, index=False)
    return path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="batos-paper-v2: run one crop end to end")
    parser.add_argument("--crop-id", required=True)
    parser.add_argument("--config", type=Path, default=Path("configs/batos_paper_v2.yaml"))
    parser.add_argument("--force", action="store_true", help="move an existing run dir aside and rerun")
    args = parser.parse_args(argv)

    t_start = time.time()
    config = load_config(args.config)
    crop = find_crop(config, args.crop_id)
    crop_dir = Path(config["paths"]["crops_root"]) / args.crop_id
    gray_path, gt_path = crop_dir / "gray.tif", crop_dir / "labels.tif"
    run_dir = Path(config["paths"]["runs_root"]) / args.crop_id
    if run_dir.exists():
        if not args.force:
            raise FileExistsError(f"{run_dir} exists; use --force to move it aside")
        backup = run_dir.with_name(f"{run_dir.name}_backup_{dt.datetime.now():%Y%m%d_%H%M%S}")
        run_dir.rename(backup)
        print(f"[run] previous run moved to {backup}")
    run_dir.mkdir(parents=True)
    timings: dict[str, float] = {}

    def log(msg):
        print(msg, flush=True)

    gray = read3d(gray_path, "gray crop")
    gt = read3d(gt_path, "GT crop")

    t = time.time()
    pre = preprocess.run(gray, config, run_dir / "01_preprocess", args.crop_id)
    timings["preprocess"] = time.time() - t
    log(f"[1/5] preprocess: R_hat={pre['scale']['estimated_grain_radius']:.2f} r={pre['openrec_radius']} "
        f"({timings['preprocess']:.0f}s)")
    openrec = read3d(pre["openrec_image"])
    selector = read3d(pre["selector_image"])

    t = time.time()
    tos = tos_select.run(selector, config, run_dir / "02_tos", args.crop_id, log=log)
    timings["tos_select"] = time.time() - t
    log(f"[2/5] tos_select: {tos['markers']['marker_count']} markers ({timings['tos_select']:.0f}s)")

    t = time.time()
    route = gate.run(Path(tos["diagnostics_csv"]), Path(tos["markers_csv"]), Path(tos["markers_tif"]),
                     openrec, config, run_dir / "03_gate", Path(pre["openrec_image"]))
    timings["gate"] = time.time() - t
    g = route["gate"]
    log(f"[3/5] gate: median={g['region_voxels_median']} q75={g['region_voxels_q75']} -> "
        f"{'TRIGGERED' if g['gate_true'] else 'not triggered'} ({route['selected_source']})")

    t = time.time()
    marker_paths = {"tos_e11_min": Path(tos["markers_tif"]), "batos_selective": Path(route["selective_markers"]),
                    "edt_otsu": Path(route["edt_otsu_markers"])}
    ws = watershed.run(openrec, selector, read3d(tos["internal_dark"]), marker_paths, config, run_dir / "04_watershed")
    timings["watershed"] = time.time() - t
    log(f"[4/5] watershed: mask fraction {ws['mask_fraction']:.3f} ({timings['watershed']:.0f}s)")

    t = time.time()
    per_crop = []
    for method, source in METHOD_MARKERS.items():
        for row in ws["rows"]:
            if row["marker_source"] != source:
                continue
            for result in evaluate.evaluate_crop(read3d(row["segmentation"]), read3d(row["normalized_marker"]),
                                                 gt, config["eval"]):
                per_crop.append({"crop_id": args.crop_id, "method": method, "size_filter": row["size_filter"],
                                 "segmentation": row["segmentation"], **result})
    eval_dir = run_dir / "05_eval"
    eval_dir.mkdir()
    per_crop_frame = pd.DataFrame(per_crop)
    per_crop_frame.to_csv(eval_dir / "per_crop_metrics.csv", index=False)
    timings["evaluate"] = time.time() - t
    timings["total"] = time.time() - t_start

    summary_cols = ["method", "size_filter", "evaluation_protocol", "n_gt_valid", "n_pred_eval", "n_pred_ignored",
                    "precision_iou050", "recall_iou050", "f1_iou025", "f1_iou050", "f1_iou075", "overseg_ratio"]
    log("[5/5] evaluation (per crop):")
    log(per_crop_frame[summary_cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    peak = peak_memory_mb()
    git = git_state()
    stamp = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    record = {
        "experiment": "batos-paper-v2", "crop_id": args.crop_id, "crop": crop,
        "command": " ".join([sys.executable, "-m", "batos_v2.run_crop", *(argv if argv is not None else sys.argv[1:])]),
        "timestamp": stamp, "git": git, "config_path": str(args.config), "config_sha256": config_hash(config),
        "config": config,
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "numpy": np.__version__, "higra": _version("higra"), "scikit-image": _version("scikit-image"),
                        "scipy": _version("scipy"), "pandas": pd.__version__, "tifffile": _version("tifffile")},
        "inputs": {"gray": str(gray_path), "gt": str(gt_path)},
        "outputs": {"run_dir": str(run_dir), "preprocess": pre, "tos_select": tos, "gate": route,
                    "watershed": {k: v for k, v in ws.items() if k != "rows"},
                    "segmentations": [r["segmentation"] for r in ws["rows"]],
                    "per_crop_metrics": str(eval_dir / "per_crop_metrics.csv")},
        "metrics": per_crop_frame.drop(columns=["segmentation"]).to_dict(orient="records"),
        "timings_seconds": timings, "peak_memory_mb": peak,
    }
    write_json(run_dir / "run.json", record)

    all_rows = []
    for r in per_crop:
        all_rows.append({
            "experiment_id": f"batos_paper_v2/{args.crop_id}/{stamp}", "method": r["method"],
            "model": "batos-paper-v2", "volume": crop["scan"], "crop": args.crop_id, "device": "cpu",
            "parameters": f"config_sha256={record['config_sha256'][:12]};size_filter={r['size_filter']};"
                          f"protocol={r['evaluation_protocol']}",
            "git_commit": git["commit"], "git_dirty": git["dirty"],
            "runtime_seconds": round(timings["total"], 1), "peak_memory_mb": round(peak, 1) if peak else None,
            "tp": r["tp_iou050"], "fp": r["fp_iou050"], "fn": r["fn_iou050"],
            "precision": r["precision_iou050"], "recall": r["recall_iou050"],
            "f1_025": r["f1_iou025"], "f1_050": r["f1_iou050"], "f1_075": r["f1_iou075"],
            "ap50": None, "ap75": None, "mean_matched_iou": r["mean_matched_iou050"],
            "n_pred": r["n_pred_eval"], "n_gt": r["n_gt_valid"], "run_json": str(run_dir / "run.json"),
        })
    log(f"[run] run.json -> {run_dir / 'run.json'}; all_runs -> {append_all_runs(all_rows, config)}")
    log(f"[run] total {timings['total']:.0f}s, peak memory {peak:.0f} MB" if peak else f"[run] total {timings['total']:.0f}s")
    return 0


def _version(dist: str) -> str | None:
    try:
        from importlib.metadata import version
        return version(dist)
    except Exception:
        return None


if __name__ == "__main__":
    raise SystemExit(main())
