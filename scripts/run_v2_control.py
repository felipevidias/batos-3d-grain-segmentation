#!/usr/bin/env python3
"""Controlled reconstruction experiment on one crop (authorised for EFRGP01_00_crop04 only).

Runs stages 2-5 of batos-paper-v2 with an explicit selector image (ToS + boundary input) and an
explicit watershed/openrec image, optionally after preprocessing overrides, and reports the
quantities the experiment log requires: pool size (MIN candidates >= tau_s), n_markers, valid
grains seeded, F1@0.50 / F1@0.75 (seed-ignore, primary) for size filter 500 and off.

Frozen and never overridable here: K, tau_s, NMS, boundary rule, gate thresholds, mask, watershed.

Usage:
  python -I scripts/run_v2_control.py --name C2_lambda15_main --crop-id EFRGP01_00_crop04 \
      --selector-image <tif> [--watershed-image <tif>] [--set openrec.radius=21 ...]
  python -I scripts/run_v2_control.py --name H1 --crop-id EFRGP01_00_crop04 --preprocess --set ...
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from batos_v2 import evaluate, gate, preprocess, tos_select, watershed
from batos_v2.common import REPO_ROOT, config_hash, load_config, read3d, write_json
from batos_v2.run_crop import git_state

FROZEN_SECTIONS = ("boundary", "markers", "gate", "watershed", "eval")
FROZEN_TOS_KEYS = ("k_total", "k_split", "min_score")
CONTROLS_ROOT = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_controls")


def apply_overrides(config: dict, items: list[str]) -> dict:
    config = copy.deepcopy(config)
    for item in items:
        key, _, raw = item.partition("=")
        section, _, field = key.partition(".")
        if section in FROZEN_SECTIONS or (section == "tos" and field in FROZEN_TOS_KEYS):
            raise ValueError(f"{key} is frozen and cannot be overridden in a control experiment")
        config.setdefault(section, {})[field] = yaml.safe_load(raw)
    return config


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--name", required=True)
    parser.add_argument("--crop-id", default="EFRGP01_00_crop04")
    parser.add_argument("--config", type=Path, default=REPO_ROOT / "configs" / "batos_paper_v2.yaml")
    parser.add_argument("--selector-image", type=Path, help="use this ToS/boundary input instead of preprocessing")
    parser.add_argument("--watershed-image", type=Path, help="default: v2 openrec image of the main run")
    parser.add_argument("--preprocess", action="store_true", help="run v2 preprocessing with --set overrides")
    parser.add_argument("--set", action="append", default=[], metavar="section.key=value")
    parser.add_argument("--hypothesis", default="", help="one-line statement recorded in result.json")
    args = parser.parse_args(argv)
    if args.crop_id != "EFRGP01_00_crop04":
        raise SystemExit("Controls are authorised for EFRGP01_00_crop04 only.")

    config = apply_overrides(load_config(args.config), args.set)
    crop_dir = Path(config["paths"]["crops_root"]) / args.crop_id
    main_run = Path(config["paths"]["runs_root"]) / args.crop_id / "01_preprocess"
    out = CONTROLS_ROOT / args.name
    if out.exists():
        raise FileExistsError(f"{out} exists; pick a new --name")
    out.mkdir(parents=True)
    t0 = time.time()
    gt = read3d(crop_dir / "labels.tif")

    if args.preprocess:
        pre = preprocess.run(read3d(crop_dir / "gray.tif"), config, out / "01_preprocess", args.crop_id)
        selector_path, ws_path = Path(pre["selector_image"]), Path(pre["openrec_image"])
    else:
        if args.selector_image is None:
            raise SystemExit("--selector-image or --preprocess is required")
        selector_path = args.selector_image
        ws_path = args.watershed_image or next(main_run.glob("*_flatten_openrec_r*_8u.tif"))
    selector, flat = read3d(selector_path), read3d(ws_path)

    tos = tos_select.run(selector, config, out / "02_tos", args.crop_id, log=lambda *a: None)
    route = gate.run(Path(tos["diagnostics_csv"]), Path(tos["markers_csv"]), Path(tos["markers_tif"]),
                     flat, config, out / "03_gate")
    ws = watershed.run(flat, selector, read3d(tos["internal_dark"]),
                       {"tos_e11_min": Path(tos["markers_tif"]), "edt_otsu": Path(route["edt_otsu_markers"])},
                       config, out / "04_watershed")

    valid, _ = evaluate.reference_sets(gt, config["eval"]["min_gt_voxels"], config["eval"]["exclude_face_touching"])
    marker_pts = np.argwhere(read3d(tos["markers_tif"]) > 0)
    seeded = len({int(v) for v in gt[tuple(marker_pts.T)] if int(v) in valid})
    rows = []
    for row in ws["rows"]:
        for r in evaluate.evaluate_crop(read3d(row["segmentation"]), read3d(row["normalized_marker"]), gt,
                                        config["eval"]):
            if r["evaluation_protocol"] == "seed_ignore":
                rows.append({"method": row["marker_source"], "size_filter": row["size_filter"],
                             "f1_iou050": r["f1_iou050"], "f1_iou075": r["f1_iou075"],
                             "precision_iou050": r["precision_iou050"], "recall_iou050": r["recall_iou050"],
                             "n_pred_eval": r["n_pred_eval"]})
    metrics = pd.DataFrame(rows)
    result = {
        "name": args.name, "hypothesis": args.hypothesis, "crop_id": args.crop_id,
        "command": " ".join([sys.executable, *sys.argv]), "git": git_state(),
        "overrides": args.set, "config_sha256": config_hash(config),
        "selector_image": str(selector_path), "watershed_image": str(ws_path),
        "pool_min_candidates": tos["tos"]["min_candidates_scored"],
        "selected_min": tos["tos"]["selected_min"], "max_min_height": tos["tos"]["max_min_height"],
        "avg_min_like_area": tos["tos"]["avg_min_like_area"], "dark_threshold": tos["boundary"]["dark_threshold"],
        "n_markers": tos["markers"]["marker_count"], "valid_gt": len(valid), "valid_grains_seeded": seeded,
        "gate": route["gate"], "metrics": metrics.to_dict(orient="records"),
        "seconds": time.time() - t0,
    }
    if args.preprocess:
        result["preprocess"] = {k: pre[k] for k in ("scale", "openrec_radius", "openrec")}
    write_json(out / "result.json", result)
    tos_rows = metrics[metrics.method == "tos_e11_min"].set_index("size_filter")
    print(json.dumps({k: result[k] for k in ("name", "pool_min_candidates", "max_min_height", "n_markers",
                                              "valid_grains_seeded", "valid_gt")}))
    print(f"gate median={route['gate']['region_voxels_median']} q75={route['gate']['region_voxels_q75']} "
          f"triggered={route['gate']['gate_true']}")
    print(metrics.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
