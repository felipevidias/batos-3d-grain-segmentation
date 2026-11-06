"""Phase 5c driver: run the extinction selector variants through the frozen downstream pipeline on one crop.

Reuses batos_v2.sensitivity.CropRun, so upstream inputs (openrec, flatten, legacy seed, background, legacy markers)
are read from the frozen run, and every downstream stage (e10 spatialise, e11 ToS markers, EDT-Otsu control,
e11 fixed-mask watershed, E23c gate, selective watershed, evaluate_v2 with markers / rule dominant / seed-ignore)
is the unmodified frozen call. Only the `true_dual` stage is replaced by batos_v2.extinction_selector.
Every stage runs as a subprocess; its peak memory is the max peak_wset over the process tree (as in frozen.py).

Outputs: C:\\BA_TOS_DATA\\runs\\extinction\\<variant>\\<crop>\\ (TIFFs kept), results/extinction/<crop>_summary.csv,
rows appended to results/all_runs.csv (one per variant x route).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
import tifffile

from .common import REPO_ROOT
from .run_crop import git_state
from .sensitivity import STAGE_SCRIPT, CropRun, RamBudget
from .extinction_selector import VARIANTS

ROOT = Path(r"C:\BA_TOS_DATA\runs\extinction")


def run_tracked(cmd: list[str], log_path: Path) -> tuple[int, float]:
    with open(log_path, "a", encoding="utf-8") as log:
        log.write(f"\n$ {subprocess.list2cmdline(cmd)}\n")
        log.flush()
        proc = psutil.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        peak = 0
        while proc.poll() is None:
            try:
                tree = [proc, *proc.children(recursive=True)]
            except psutil.Error:
                tree = [proc]
            for p in tree:
                try:
                    info = p.memory_info()
                    peak = max(peak, getattr(info, "peak_wset", 0) or info.rss)
                except psutil.Error:
                    pass
            time.sleep(0.2)
    return proc.returncode, peak / 2 ** 20


class ExtinctionRun(CropRun):
    def run(self, name: str, stage: str, args: list) -> None:
        if stage == "true_dual":
            image, out_dir = args[0], args[args.index("--out-dir") + 1]
            cmd = [sys.executable, "-I", "-m", "batos_v2.extinction_selector", str(image), "--out-dir", str(out_dir),
                   "--variant", str(self.value), "--k", "3600"]
            script = "batos_v2.extinction_selector"
        else:
            script = str(self.scripts[STAGE_SCRIPT[stage]])
            cmd = [sys.executable, "-I", script, *[str(a) for a in args]]
        t0 = time.time()
        rc, peak_mb = run_tracked(cmd, self.log)
        self.steps.append({"step": name, "script": script, "seconds": round(time.time() - t0, 1), "rc": rc,
                           "peak_memory_mb": round(peak_mb, 1)})
        print(f"[ext] {self.value} {name}: rc={rc} {time.time() - t0:.0f}s peak {peak_mb:.0f} MB", flush=True)
        if rc != 0:
            raise RuntimeError(f"{name} failed rc={rc}; see {self.log}")


def gt_path(crop: str) -> Path:
    meta = json.loads((REPO_ROOT / "data/metadata/crops_sibgrapi.json").read_text(encoding="utf-8"))
    for c in meta["crops"]:
        if c["crop_id"] == crop:
            return Path(c["gt_path"])
    raise KeyError(crop)


def valid_gt(gt: np.ndarray) -> np.ndarray:
    """Valid reference labels: >= 1000 voxels and not touching a face (evaluate_v2 defaults)."""
    labels, counts = np.unique(gt[gt > 0], return_counts=True)
    face = set()
    for s in (gt[0], gt[-1], gt[:, 0], gt[:, -1], gt[:, :, 0], gt[:, :, -1]):
        face |= set(np.unique(s).tolist())
    return np.array([l for l, c in zip(labels, counts) if c >= 1000 and l not in face])


def seeded(markers: np.ndarray, gt: np.ndarray, valid: np.ndarray) -> tuple[int, int]:
    """(valid GT grains containing >= 1 marker voxel, number of marker labels). GT used for diagnosis only."""
    hit = np.unique(gt[(markers > 0) & (gt > 0)])
    return int(np.isin(valid, hit).sum()), int(len(np.unique(markers[markers > 0])))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Phase 5c extinction selector on one crop")
    ap.add_argument("--crop", default="EFRGP01_00_crop04")
    ap.add_argument("--variants", nargs="*", default=list(VARIANTS))
    ap.add_argument("--ram-gb", type=float, default=6.0)
    args = ap.parse_args(argv)
    git = git_state()
    rows, all_runs = [], []
    for v in args.variants:
        out = ROOT / v / args.crop
        param = {"name": f"extinction_{v}", "group": "extinction", "kind": "flag", "stage": "true_dual",
                 "flag": "--variant", "frozen": "weighted_true_dual", "type": "str"}
        t0 = time.time()
        if not (out / "DONE").exists():
            if out.exists():
                raise RuntimeError(f"partial run folder exists, inspect and remove it first: {out}")
            run = ExtinctionRun(args.crop, param, "replace_selector", v, out, RamBudget(args.ram_gb), keep_tifs=True)
            run.note = "true_dual replaced by batos_v2.extinction_selector (flags['true_dual'] are NOT used)"
            run.execute()
        rec = json.loads((out / "run.json").read_text(encoding="utf-8"))
        rec["git"] = git
        (out / "run.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
        info = json.loads((out / "e11" / "selector" / "extinction_selector_info.json").read_text(encoding="utf-8"))
        metrics = pd.read_csv(out / "metrics.csv")
        gt = tifffile.imread(gt_path(args.crop))
        valid = valid_gt(gt)
        peak = max(s.get("peak_memory_mb", 0) for s in rec["steps"])
        secs = sum(s["seconds"] for s in rec["steps"])
        for route, mpath in (("BA-TOS_ToS", out / "e11" / "watershed" / "normalized_markers" / "tos_e11_min.tif"),
                             ("BA-TOS_sel", out / "e23c" / "watershed" / "normalized_markers" / "batos_selective.tif")):
            m = metrics[metrics["route"] == route].iloc[0]
            n_seed, n_mark = seeded(tifffile.imread(mpath), gt, valid)
            rows.append({"crop_id": args.crop, "variant": v, "route": route,
                         "dark_extrema": info["num_dark_extrema"],
                         "dark_extrema_touching_face": info["num_dark_extrema_touching_face"],
                         "pool_size": info["pool_size"], "selected_K": info["num_selected"],
                         "rep_area_median": info["representative_area_median"],
                         "rep_area_q75": info["representative_area_q75"],
                         "n_markers": n_mark, "valid_gt": len(valid), "grains_seeded": n_seed,
                         "seeded_fraction": n_seed / max(1, len(valid)),
                         "gate_triggered": bool(m["gate_triggered"]), "gate_mv": m["region_voxels_median"],
                         "gate_qv": m["region_voxels_q75"], "n_pred_evaluable": m["prediction_count_evaluable"],
                         "precision_050": m["precision_iou050"], "recall_050": m["recall_iou050"],
                         "f1_025": m["f1_iou025"], "f1_050": m["f1_iou050"], "f1_075": m["f1_iou075"],
                         "selector_seconds": info["seconds"], "pipeline_seconds": secs,
                         "max_stage_peak_mb": peak})
            all_runs.append({
                "experiment_id": f"phase5c_extinction_{v}/{args.crop}", "method": route,
                "model": f"batos-paper-v2 with extinction selector ({v})", "volume": args.crop[:10],
                "crop": args.crop, "device": "cpu (AMD Ryzen 7 5700X)",
                "parameters": f"selector=extinction_{v};K=3600;size_filter=none;protocol=seed_ignore;NMS=6",
                "runtime_seconds": secs, "peak_memory_mb": peak, "tp": m["tp_iou050"], "fp": m["fp_iou050"],
                "fn": m["fn_iou050"], "precision": m["precision_iou050"], "recall": m["recall_iou050"],
                "f1_025": m["f1_iou025"], "f1_050": m["f1_iou050"], "f1_075": m["f1_iou075"], "ap50": None,
                "ap75": None, "mean_matched_iou": m.get("mean_matched_iou050"),
                "n_pred": m["prediction_count_evaluable"], "n_gt": m["valid_gt_count"],
                "git_commit": git["commit"], "run_json": str(out / "run.json"), "memory_run_json": str(out / "run.json")})
        print(f"[ext] {v} done in {time.time() - t0:.0f}s", flush=True)
    res = REPO_ROOT / "results" / "extinction"
    res.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(res / f"{args.crop}_summary.csv", index=False)
    path = REPO_ROOT / "results" / "all_runs.csv"
    frame = pd.read_csv(path)
    frame = frame[~frame["experiment_id"].astype(str).str.startswith("phase5c_extinction_")
                  | ~frame["crop"].eq(args.crop)]
    pd.concat([frame, pd.DataFrame(all_runs)], ignore_index=True).to_csv(path, index=False)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)
    print(f"[ext] driver peak {psutil.Process().memory_info().peak_wset / 2**20:.0f} MB", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
