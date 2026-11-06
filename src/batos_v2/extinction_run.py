"""Phase 5c driver: run extinction selector variants through the frozen downstream pipeline.

Reuses batos_v2.sensitivity.CropRun, so upstream inputs (openrec, flatten, legacy seed, background, legacy markers)
are read from the frozen run, and every downstream stage (e10 spatialise, e11 ToS markers, EDT-Otsu control,
e11 fixed-mask watershed, E23c gate, selective watershed, evaluate_v2 with markers / rule dominant / seed-ignore)
is the unmodified frozen call. Only the `true_dual` stage is replaced by batos_v2.extinction_selector.
Every stage runs as a subprocess; its peak memory is the max peak_wset over the process tree (as in frozen.py).

volume_unique* variants get --internal-source = the frozen internal-dark mask of the crop; volume_unique_scale gets
--scale-radius = the frozen e13 r_est (tag "r_est") or, with --scale-json, the granulometric estimate (tag
"r_<key>").

Outputs: C:\\BA_TOS_DATA\\runs\\extinction\\<variant>[_<tag>]\\<crop>\\ (TIFFs kept),
results/extinction/summary.csv (one row per crop x variant x route, upserted), rows upserted in results/all_runs.csv.
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
from .sensitivity import FROZEN, STAGE_SCRIPT, CropRun, RamBudget
from .extinction_selector import VARIANTS

ROOT = Path(r"C:\BA_TOS_DATA\runs\extinction")
RESULTS = REPO_ROOT / "results" / "extinction"


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


class TrackedRun(CropRun):
    """CropRun whose stages are tracked for peak memory; `replace` maps stage -> command builder."""
    replace: dict = {}

    def run(self, name: str, stage: str, args: list) -> None:
        if stage in self.replace:
            cmd, script = self.replace[stage](self, args)
        else:
            script = str(self.scripts[STAGE_SCRIPT[stage]])
            cmd = [sys.executable, "-I", script, *[str(a) for a in args]]
        t0 = time.time()
        rc, peak_mb = run_tracked(cmd, self.log)
        self.steps.append({"step": name, "script": script, "seconds": round(time.time() - t0, 1), "rc": rc,
                           "peak_memory_mb": round(peak_mb, 1)})
        print(f"[run] {self.crop} {self.value} {name}: rc={rc} {time.time() - t0:.0f}s peak {peak_mb:.0f} MB",
              flush=True)
        if rc != 0:
            raise RuntimeError(f"{name} failed rc={rc}; see {self.log}")


class ExtinctionRun(TrackedRun):
    extra_args: list = []

    def _selector(self, args):
        image, out_dir = args[0], args[args.index("--out-dir") + 1]
        cmd = [sys.executable, "-I", "-m", "batos_v2.extinction_selector", str(image), "--out-dir", str(out_dir),
               "--variant", str(self.value), "--k", "3600", *map(str, self.extra_args)]
        return cmd, "batos_v2.extinction_selector"

    replace = {"true_dual": _selector}


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


def marker_diagnosis(markers: np.ndarray, gt: np.ndarray, valid: np.ndarray) -> dict:
    """Markers by dominant GT label: valid grains seeded / split (>= 2 markers), markers in background / excluded."""
    labs = np.unique(markers[markers > 0])
    dom = np.zeros(len(labs), dtype=np.int64)
    if len(labs):
        idx = {int(l): i for i, l in enumerate(labs)}
        mk = markers[markers > 0]
        g = gt[markers > 0].astype(np.int64)
        frame = pd.DataFrame({"m": mk, "g": g})
        frame = frame[frame["g"] > 0]
        if len(frame):
            top = frame.groupby(["m", "g"]).size().reset_index(name="n").sort_values(["m", "n", "g"],
                                                                                     ascending=[True, False, True])
            top = top.drop_duplicates("m")
            for m, gl in zip(top["m"], top["g"]):
                dom[idx[int(m)]] = gl
    in_valid = np.isin(dom, valid)
    per = pd.Series(dom[in_valid]).value_counts()
    return {"n_markers": int(len(labs)), "markers_in_valid_grains": int(in_valid.sum()),
            "markers_in_background": int((dom == 0).sum()),
            "markers_in_excluded_grains": int(((dom > 0) & ~in_valid).sum()),
            "grains_seeded": int(len(per)), "grains_split": int((per >= 2).sum()),
            "extra_markers_from_splits": int((per - 1).sum())}


def summarize(out: Path, crop: str, label: str, info: dict | None, git: dict) -> tuple[list[dict], list[dict]]:
    rec = json.loads((out / "run.json").read_text(encoding="utf-8"))
    metrics = pd.read_csv(out / "metrics.csv")
    gt = tifffile.imread(gt_path(crop))
    valid = valid_gt(gt)
    peak = max(s.get("peak_memory_mb", 0) for s in rec["steps"])
    secs = sum(s["seconds"] for s in rec["steps"])
    rows, runs = [], []
    for route, mpath in (("BA-TOS_ToS", out / "e11" / "watershed" / "normalized_markers" / "tos_e11_min.tif"),
                         ("BA-TOS_sel", out / "e23c" / "watershed" / "normalized_markers" / "batos_selective.tif")):
        m = metrics[metrics["route"] == route].iloc[0]
        diag = marker_diagnosis(tifffile.imread(mpath), gt, valid)
        row = {"crop_id": crop, "variant": label, "route": route, "valid_gt": len(valid), **diag,
               "gate_triggered": bool(m["gate_triggered"]), "gate_mv": m["region_voxels_median"],
               "gate_qv": m["region_voxels_q75"], "n_pred_evaluable": m["prediction_count_evaluable"],
               "precision_050": m["precision_iou050"], "recall_050": m["recall_iou050"], "f1_025": m["f1_iou025"],
               "f1_050": m["f1_iou050"], "f1_075": m["f1_iou075"], "openrec_radius": rec.get("openrec_radius"),
               "pipeline_seconds": secs, "max_stage_peak_mb": peak}
        if info:
            row |= {k: info.get(k) for k in ("num_dark_extrema", "pool_before_rules", "pool_in_scale_window",
                                              "internal_components", "pool_without_component",
                                              "components_with_minimum", "pool_size", "num_selected",
                                              "representative_area_median", "representative_area_q75",
                                              "scale_radius")}
        rows.append(row)
        runs.append({
            "experiment_id": f"phase5_{label}/{crop}", "method": route, "model": f"batos-paper-v2 ({label})",
            "volume": crop[:10], "crop": crop, "device": "cpu (AMD Ryzen 7 5700X)",
            "parameters": f"{label};size_filter=none;protocol=seed_ignore;NMS=6;openrec_r={rec.get('openrec_radius')}",
            "runtime_seconds": secs, "peak_memory_mb": peak, "tp": m["tp_iou050"], "fp": m["fp_iou050"],
            "fn": m["fn_iou050"], "precision": m["precision_iou050"], "recall": m["recall_iou050"],
            "f1_025": m["f1_iou025"], "f1_050": m["f1_iou050"], "f1_075": m["f1_iou075"], "ap50": None, "ap75": None,
            "mean_matched_iou": m.get("mean_matched_iou050"), "n_pred": m["prediction_count_evaluable"],
            "n_gt": m["valid_gt_count"], "git_commit": git["commit"], "run_json": str(out / "run.json"),
            "memory_run_json": str(out / "run.json")})
    return rows, runs


def upsert(path: Path, rows: list[dict], keys: list[str]) -> None:
    new = pd.DataFrame(rows)
    if path.exists():
        old = pd.read_csv(path)
        if all(k in old.columns for k in keys):
            merged = old.merge(new[keys].drop_duplicates(), on=keys, how="left", indicator=True)
            old = old[(merged["_merge"] == "left_only").to_numpy()]
        new = pd.concat([old, new], ignore_index=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    new.to_csv(path, index=False)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Phase 5c extinction selector on given crops")
    ap.add_argument("--crops", nargs="*", default=["EFRGP01_00_crop04"])
    ap.add_argument("--variants", nargs="*", default=["volume", "dynamics", "product"], choices=VARIANTS)
    ap.add_argument("--scale-json", type=Path, help="results/scale/{crop}_granulo.json pattern dir for r")
    ap.add_argument("--scale-key", default="r_ps_peak")
    ap.add_argument("--scale-on", choices=("region_voxels", "volume_extinction"), default="region_voxels")
    ap.add_argument("--ram-gb", type=float, default=6.0)
    args = ap.parse_args(argv)
    git = git_state()
    rows, runs = [], []
    for crop in args.crops:
        auto = json.loads((FROZEN / crop / "frozen_run.json").read_text(encoding="utf-8"))["auto_params"]
        for v in args.variants:
            extra, label = [], f"extinction_{v}"
            if v.startswith("volume_unique"):
                internal = json.loads((FROZEN / crop / "frozen_run.json").read_text(encoding="utf-8"))["outputs"][
                    "internal"]
                extra += ["--internal-source", internal]
            if v == "volume_unique_scale":
                if args.scale_json:
                    r = json.loads((args.scale_json / f"{crop}_granulo.json").read_text(encoding="utf-8"))[
                        args.scale_key]
                    tag = args.scale_key
                else:
                    r, tag = float(auto["estimated_grain_radius"]), "r_est"
                extra += ["--scale-radius", r, "--scale-on", args.scale_on]
                label = f"extinction_{v}_{tag}" + ("_greyvol" if args.scale_on == "volume_extinction" else "")
            out = ROOT / label.removeprefix("extinction_") / crop
            if not (out / "DONE").exists():
                if out.exists():
                    raise RuntimeError(f"partial run folder exists, inspect and remove it first: {out}")
                param = {"name": label, "group": "extinction", "kind": "flag", "stage": "true_dual",
                         "flag": "--variant", "frozen": "weighted_true_dual", "type": "str"}
                run = ExtinctionRun(crop, param, "replace_selector", v, out, RamBudget(args.ram_gb), keep_tifs=True)
                run.extra_args = extra
                run.note = f"true_dual replaced by batos_v2.extinction_selector {v} {extra} (flags['true_dual'] unused)"
                run.execute()
                rec = json.loads((out / "run.json").read_text(encoding="utf-8"))
                rec["git"] = git
                (out / "run.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
            info = json.loads((out / "e11" / "selector" / "extinction_selector_info.json").read_text(encoding="utf-8"))
            r1, r2 = summarize(out, crop, label, info, git)
            rows += r1
            runs += r2
    upsert(RESULTS / "summary.csv", rows, ["crop_id", "variant", "route"])
    upsert(REPO_ROOT / "results" / "all_runs.csv", runs, ["experiment_id", "method"])
    print(pd.DataFrame(rows)[["crop_id", "variant", "route", "n_markers", "grains_seeded", "grains_split",
                              "gate_triggered", "gate_mv", "gate_qv", "f1_050", "f1_075",
                              "max_stage_peak_mb"]].to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
