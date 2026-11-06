"""Phase 5a: global one-factor perturbation sensitivity of the frozen BA-TOS pipeline on all 12 crops.

For each parameter of configs/sensitivity_global.yaml and each perturbed value, the frozen pipeline is
re-run on the 12 crops FROM THE FIRST STAGE THAT USES THE PARAMETER. Upstream stage outputs are read
(never written) from the frozen runs in C:\\BA_TOS_DATA\\runs\\batos_paper_v2_frozen. Every stage is the
unmodified recovered script called as a subprocess (python -I), with the frozen flags of batos_v2.frozen,
except for the one perturbation:
  kind flag     : one CLI flag of the stage's command line is set to the perturbed value (appended if the
                  frozen call relies on the script's default);
  kind literal  : one hard-coded constant: the script is COPIED to <run>/_patched/ with exactly one source
                  substitution (old -> new; fails loudly unless `old` occurs exactly once). The frozen file
                  is never touched; the patched copy's SHA-256 and the substitution are recorded;
  kind scale    : openrec radius = max(floor, round(factor * R_hat)), R_hat from the frozen GT-free
                  e13 scale estimate (frozen factor 1.5, floor 8; identity reproduces the frozen radius);
  kind eval     : evaluation protocol only (evaluate_v2 kwargs) on the frozen segmentations.
Routes evaluated: BA-TOS_ToS (tos_e11_min, e11 watershed) and BA-TOS_sel (E23c gate route), with
batos_v2.evaluate_v2 using the watershed's normalized markers (rule dominant; == e09b seed_ignore, Phase 1d).
Output per run: C:\\BA_TOS_DATA\\runs\\sensitivity\\<param>_<value>\\<crop>\\{metrics.csv, run.json, commands.log};
large TIFF intermediates are deleted after evaluation (SHA-256 of both segmentations kept in run.json).
Parallel: crops (and perturbations) run in a thread pool; every subprocess acquires a RAM budget weighted
by the stage's measured peak (frozen_run.json of the memory runs), total budget --ram-gb (default 22).
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import yaml

from .common import REPO_ROOT
from .evaluate_v2 import evaluate
from .frozen import BACKGROUND, MARKER_FLAGS, MARKERS_LEGACY, SCRIPTS, TOS_SEED_LEGACY, TRUE_DUAL
from .run_crop import git_state

FROZEN = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen")
DEFAULT_ROOT = Path(r"C:\BA_TOS_DATA\runs\sensitivity")
CROPS = [f"EFRGP01_0{s}_crop0{c}" for s in (0, 1) for c in (2, 3, 4, 5, 6, 7)]
STAGES = ["openrec", "flatten", "tos_seed", "background", "markers_legacy", "true_dual", "spatialize",
          "markers_tos", "edt_otsu", "ws_tos", "route", "ws_sel", "evaluate"]
DEPS = {"openrec": [], "flatten": ["openrec"], "tos_seed": ["flatten"],
        "background": ["openrec", "flatten", "tos_seed"], "markers_legacy": ["openrec", "flatten", "background"],
        "true_dual": ["flatten"], "spatialize": ["flatten", "true_dual", "background"],
        "markers_tos": ["spatialize", "background"], "edt_otsu": ["openrec", "markers_tos"],
        "ws_tos": ["openrec", "markers_tos", "markers_legacy"], "route": ["spatialize", "markers_tos", "edt_otsu"],
        "ws_sel": ["openrec", "markers_tos", "ws_tos", "route"], "evaluate": ["ws_tos", "ws_sel"]}
STAGE_SCRIPT = {"openrec": "openrec_v2", "flatten": "batos_preprocess", "tos_seed": "batos_tos_seed",
                "background": "batos_background", "markers_legacy": "batos_markers", "true_dual": "true_dual",
                "spatialize": "spatialize", "markers_tos": "e11_markers", "edt_otsu": "e11_controls",
                "ws_tos": "e11_watershed", "route": "e23c_route", "ws_sel": "e11_watershed"}
STAGE_RAM_GB = {"openrec": 4.8, "tos_seed": 3.3, "true_dual": 3.3, "spatialize": 3.3}  # others 0.7
FLATTEN = ["--lambda-xy", "1"]
EDT_OTSU = ["--mode", "otsu_edt", "--nms-radius", "6"]
WS_TOS_POLICY = ["--outside-policy", "drop"]
BASE_FLAGS = {"tos_seed": TOS_SEED_LEGACY, "background": BACKGROUND, "markers_legacy": MARKERS_LEGACY,
              "true_dual": TRUE_DUAL, "markers_tos": ["--polarity", "min", *MARKER_FLAGS], "flatten": FLATTEN,
              "edt_otsu": EDT_OTSU, "ws_tos": WS_TOS_POLICY, "ws_sel": WS_TOS_POLICY}


def downstream(stage: str) -> set[str]:
    dirty = {stage}
    changed = True
    while changed:
        changed = False
        for s, deps in DEPS.items():
            if s not in dirty and dirty & set(deps):
                dirty.add(s)
                changed = True
    return dirty


def set_flag(flags: list[str], flag: str, value) -> list[str]:
    out = [str(f) for f in flags]
    if flag in out:
        i = out.index(flag)
        if value is True or value is False:
            raise ValueError(f"boolean flag {flag} given with a value in the frozen call")
        out[i + 1] = str(value)
    elif value is True:
        out.append(flag)
    elif value is False:
        pass
    else:
        out += [flag, str(value)]
    return out


def unset_store_true(flags: list[str], flag: str) -> list[str]:
    return [f for f in flags if f != flag]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha256_array(path: Path) -> str:
    return hashlib.sha256(np.ascontiguousarray(tifffile.imread(path)).tobytes()).hexdigest()


def fmt_value(v) -> str:
    if isinstance(v, bool):
        return str(v).lower()
    if isinstance(v, float):
        s = f"{v:.6g}"
    else:
        s = str(v)
    return s.replace("-", "m").replace(".", "p")


# --------------------------------------------------------------------------- RAM budget
class RamBudget:
    def __init__(self, gb: float):
        self.total, self.used = gb, 0.0
        self.cv = threading.Condition()

    def acquire(self, gb: float):
        with self.cv:
            while self.used + gb > self.total + 1e-9:
                self.cv.wait()
            self.used += gb

    def release(self, gb: float):
        with self.cv:
            self.used -= gb
            self.cv.notify_all()


class PeakMonitor(threading.Thread):
    """System-wide used RAM, sampled every second (includes everything else running on the machine)."""

    def __init__(self):
        super().__init__(daemon=True)
        import psutil
        self.psutil = psutil
        self.base = psutil.virtual_memory().used
        self.peak = self.base
        self.stop = threading.Event()

    def run(self):
        while not self.stop.is_set():
            self.peak = max(self.peak, self.psutil.virtual_memory().used)
            time.sleep(1.0)


# --------------------------------------------------------------------------- perturbations
def perturbed_values(p: dict) -> list[tuple[str, object]]:
    """(direction, value) pairs. Rules: float -> x0.5, x2; int >= 4 -> round(x0.5), round(x2);
    int < 4 -> -1, +1; bool -> flip; explicit `values` in the config override; domain clipped via lo/hi."""
    v, typ = p["frozen"], p["type"]
    if "values" in p:
        pairs = [(d, val) for d, val in p["values"].items()]
    elif typ == "bool":
        pairs = [("flip", not v)]
    elif typ == "int" and abs(v) >= 4:
        pairs = [("x0.5", int(round(v * 0.5))), ("x2", int(round(v * 2)))]
    elif typ == "int":
        pairs = [("-1", v - 1), ("+1", v + 1)]
    else:
        pairs = [("x0.5", v * 0.5), ("x2", v * 2.0)]
    out = []
    for d, val in pairs:
        if typ != "bool" and val is not None:
            if "lo" in p:
                val = max(p["lo"], val)
            if "hi" in p:
                val = min(p["hi"], val)
            if typ == "int":
                val = int(val)
        if val is None or val == v:
            continue
        out.append((d, val))
    return out


def load_params(config: Path) -> list[dict]:
    cfg = yaml.safe_load(Path(config).read_text(encoding="utf-8"))
    params = []
    for group, items in cfg["groups"].items():
        for p in items:
            params.append({"group": group, **p})
    names = [p["name"] for p in params]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        raise ValueError(f"duplicate parameter names: {dup}")
    return params


# --------------------------------------------------------------------------- one crop, one perturbation
class CropRun:
    def __init__(self, crop: str, param: dict | None, direction: str, value, out: Path, budget: RamBudget,
                 keep_tifs: bool = False):
        self.crop, self.param, self.direction, self.value = crop, param, direction, value
        self.out, self.budget, self.keep_tifs = out, budget, keep_tifs
        self.log = out / "commands.log"
        self.steps: list[dict] = []
        self.scripts = dict(SCRIPTS)
        self.patch_record = None

    # frozen paths
    def frozen_paths(self) -> dict[str, Path]:
        f = FROZEN / self.crop
        rec = json.loads((f / "frozen_run.json").read_text(encoding="utf-8"))
        o = rec["outputs"]
        pipe = f / "batos_openrec_auto_minparam"
        seeds = sorted((pipe / "02_tos_seed").glob("*_tosconop_autodual_seed_raw.tif"))
        paths = {"gray": Path(rec["inputs"]["gray"]), "gt": Path(rec["inputs"]["gt"]),
                 "flat": Path(o["flat_image"]), "selector": Path(o["selector"]),
                 "raw_seed": seeds[0] if len(seeds) == 1 else None, "internal": Path(o["internal"]),
                 "baseline": Path(o["baseline_markers"]), "cand_csv": Path(o["candidates"]),
                 "diag": Path(o["diagnostics"]),
                 "tos_dir": f / "e11" / "markers_tos_min", "edt_dir": f / "e11" / "markers_edt_otsu",
                 "ws_tos": f / "e11" / "watershed", "route_dir": f / "e23c" / "route", "ws_sel": f / "e23c" / "watershed"}
        for k, v in paths.items():
            if v is None or not Path(v).exists():
                raise FileNotFoundError(f"{self.crop} frozen {k}: {v}")
        self.auto = rec["auto_params"]
        self.frozen_radius = int(rec["openrec_radius"])
        return paths

    def run(self, name: str, stage: str, args: list) -> None:
        script = self.scripts[STAGE_SCRIPT[stage]]
        cmd = [sys.executable, "-I", str(script), *[str(a) for a in args]]
        gb = STAGE_RAM_GB.get(stage, 0.7)
        self.budget.acquire(gb)
        t0 = time.time()
        try:
            with open(self.log, "a", encoding="utf-8") as log:
                log.write(f"\n$ {subprocess.list2cmdline(cmd)}\n")
                log.flush()
                rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
        finally:
            self.budget.release(gb)
        self.steps.append({"step": name, "script": str(script), "seconds": round(time.time() - t0, 1), "rc": rc})
        if rc != 0:
            raise RuntimeError(f"{name} failed rc={rc}; see {self.log}")

    def patch_literal(self, stage: str, olds: list[str], news: list[str], count: int = 1) -> None:
        key = STAGE_SCRIPT[stage]
        src = Path(SCRIPTS[key])
        text = src.read_text(encoding="utf-8")
        for old, new in zip(olds, news, strict=True):
            n = text.count(old)
            if n != count:
                raise RuntimeError(f"literal patch: {old!r} occurs {n} times in {src} (expected {count})")
            text = text.replace(old, new)
        patched_dir = self.out / "_patched"
        patched_dir.mkdir(parents=True, exist_ok=True)
        dst = patched_dir / src.name
        dst.write_text(text, encoding="utf-8")
        self.scripts[key] = dst
        self.patch_record = {"source": str(src), "source_sha256": sha256_file(src), "patched": str(dst),
                             "patched_sha256": sha256_file(dst), "old": olds, "new": news, "count": count}

    def execute(self) -> dict:
        self.out.mkdir(parents=True, exist_ok=True)
        fz = self.frozen_paths()
        p, v = self.param, self.value
        flags = {k: list(map(str, x)) for k, x in BASE_FLAGS.items()}
        radius = self.frozen_radius
        eval_kwargs = {}
        dirty: set[str] = set()
        if p is not None:
            kind, stage = p["kind"], p.get("stage")
            if kind == "flag":
                if p["type"] == "bool" and v is False:
                    flags[stage] = unset_store_true(flags.get(stage, []), p["flag"])
                else:
                    flags[stage] = set_flag(flags.get(stage, []), p["flag"], v)
                dirty = downstream(stage)
            elif kind == "literal":
                olds = p["old"] if isinstance(p["old"], list) else [p["old"]]
                news = p["new"] if isinstance(p["new"], list) else [p["new"]]
                comp = 100 - v if isinstance(v, (int, float)) and not isinstance(v, bool) else None
                news = [n.format(v=v, c=comp) for n in news]
                self.patch_literal(stage, olds, news, int(p.get("count", 1)))
                dirty = downstream(stage)
            elif kind == "scale":
                factor = v if p["name"] == "openrec_radius_factor" else 1.5
                floor = v if p["name"] == "openrec_radius_floor" else 8
                radius = max(int(floor), int(round(factor * float(self.auto["estimated_grain_radius"]))))
                dirty = downstream("openrec") if radius != self.frozen_radius else {"evaluate"}
            elif kind == "eval":
                eval_kwargs = {p["kwarg"]: v}
                dirty = {"evaluate"}
            else:
                raise ValueError(kind)
        else:
            dirty = set(STAGES)  # identity run: everything from scratch
        if p is not None and p["kind"] == "scale" and "openrec" not in dirty:
            self.note = f"radius unchanged ({radius}); frozen segmentations re-evaluated"
        cur = dict(fz)
        o = self.out
        cid = self.crop
        if "openrec" in dirty:
            prep = o / "01_openrec"
            self.run("openrec", "openrec", [cur["gray"], "--out-dir", prep, "--prefix", cid, "--radius", radius])
            cur["flat"] = one(prep, "*_flatten_openrec_r*_8u.tif")
        pipe = o / "pipe"
        stem = cur["flat"].stem
        if "flatten" in dirty:
            self.run("flatten", "flatten", [cur["flat"], "--out-dir", pipe / "01_preprocess", "--prefix", stem,
                                            *flags["flatten"]])
            cur["selector"] = one(pipe / "01_preprocess", f"{stem}_flatten_lambda*_8u.tif")
        if "tos_seed" in dirty:
            self.run("tos_seed", "tos_seed", [cur["selector"], "--out-dir", pipe / "02_tos_seed", *flags["tos_seed"]])
            cur["raw_seed"] = one(pipe / "02_tos_seed", "*_tosconop_autodual_seed_raw.tif")
        if "background" in dirty:
            self.run("background", "background", [cur["flat"], "--flattened", cur["selector"], "--seed", cur["raw_seed"],
                                                  "--out-dir", pipe / "03_background", "--prefix", stem,
                                                  *flags["background"]])
            cur["internal"] = one(pipe / "03_background", "*_internal_dark_candidates.tif")
        if "markers_legacy" in dirty:
            self.run("markers_legacy", "markers_legacy",
                     [cur["flat"], "--flattened", cur["selector"], "--internal-dark", cur["internal"],
                      "--out-dir", pipe / "04_markers", "--prefix", stem, *flags["markers_legacy"]])
            cur["baseline"] = one(pipe / "04_markers", "*_batos_marker_components_raw.tif")
        e11 = o / "e11"
        if "true_dual" in dirty:
            self.run("true_dual", "true_dual", [cur["selector"], "--out-dir", e11 / "selector", *flags["true_dual"]])
            cur["cand_csv"] = one(e11 / "selector", "*_tosconop_autodual_candidates.csv")
        if "spatialize" in dirty:
            self.run("spatialize", "spatialize", ["--image", cur["selector"], "--candidate-csv", cur["cand_csv"],
                                                  "--internal-source", cur["internal"], "--out-dir", e11 / "spatial"])
            cur["diag"] = e11 / "spatial" / "candidate_diagnostics_no_gt.csv"
        if "markers_tos" in dirty:
            self.run("markers_tos", "markers_tos", ["--diagnostic", cur["diag"], "--internal-source", cur["internal"],
                                                    *flags["markers_tos"], "--out-dir", e11 / "markers_tos_min"])
            cur["tos_dir"] = e11 / "markers_tos_min"
        tos_marker = cur["tos_dir"] / "markers_uint16.tif"
        if "edt_otsu" in dirty:
            self.run("edt_otsu", "edt_otsu", ["--image", cur["flat"], *flags["edt_otsu"],
                                              "--match-count-from", tos_marker, "--out-dir", e11 / "markers_edt_otsu"])
            cur["edt_dir"] = e11 / "markers_edt_otsu"
        if "ws_tos" in dirty:
            self.run("ws_tos", "ws_tos", ["--image", cur["flat"], "--mask-marker", tos_marker,
                                          "--mask-marker", cur["baseline"], *flags["ws_tos"],
                                          "--marker", f"tos_e11_min={tos_marker}",
                                          "--marker", f"baseline_existing={cur['baseline']}",
                                          "--out-dir", e11 / "watershed"])
            cur["ws_tos"] = e11 / "watershed"
        e23 = o / "e23c"
        if "route" in dirty:
            self.run("route", "route", ["--candidate-diagnostics", cur["diag"],
                                        "--tos-marker-csv", cur["tos_dir"] / "markers.csv", "--tos-marker-tiff", tos_marker,
                                        "--edt-marker-tiff", cur["edt_dir"] / "markers_uint16.tif",
                                        "--crop-id", cid, "--out-dir", e23 / "route"])
            cur["route_dir"] = e23 / "route"
        if "ws_sel" in dirty:
            self.run("ws_sel", "ws_sel", ["--image", cur["flat"], "--mask-marker", tos_marker,
                                          "--mask-marker", cur["ws_tos"] / "normalized_markers" / "baseline_existing.tif",
                                          *flags["ws_sel"],
                                          "--marker", f"batos_selective={cur['route_dir'] / 'batos_selective_markers_uint16.tif'}",
                                          "--out-dir", e23 / "watershed"])
            cur["ws_sel"] = e23 / "watershed"

        # evaluate (always)
        gt = tifffile.imread(cur["gt"])
        rows, shas = [], {}
        gate = json.loads((cur["route_dir"] / "route_manifest.json").read_text(encoding="utf-8"))
        gate = gate.get("gate", gate)
        for route, ws, label in (("BA-TOS_ToS", cur["ws_tos"], "tos_e11_min"),
                                 ("BA-TOS_sel", cur["ws_sel"], "batos_selective")):
            seg_path = ws / "segmentations" / f"{label}_distance_c1.tif"
            shas[route] = sha256_array(seg_path)
            m, _, _ = evaluate(tifffile.imread(seg_path), gt,
                               tifffile.imread(ws / "normalized_markers" / f"{label}.tif"),
                               rule="dominant", morph_iou=None, **eval_kwargs)
            rows.append({"crop_id": self.crop, "route": route, **m})
        gate_info = {"gate_triggered": bool(gate.get("gate_triggered")),
                     "region_voxels_median": gate.get("region_voxels_median"),
                     "region_voxels_q75": gate.get("region_voxels_q75")}
        pd.DataFrame(rows).assign(**gate_info).to_csv(o / "metrics.csv", index=False)
        record = {"crop_id": self.crop, "param": None if p is None else p["name"], "direction": self.direction,
                  "value": v, "frozen": None if p is None else p["frozen"], "dirty_stages": sorted(dirty, key=STAGES.index),
                  "openrec_radius": radius, "flags": flags, "eval_kwargs": eval_kwargs, "patch": self.patch_record,
                  "segmentation_sha256": shas, "gate": gate_info, "steps": self.steps,
                  "timestamp": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
                  "note": getattr(self, "note", "")}
        (o / "run.json").write_text(json.dumps(record, indent=2, default=str), encoding="utf-8")
        if not self.keep_tifs:
            for t in list(o.rglob("*.tif")) + list(o.rglob("*.npy")):
                t.unlink()
        (o / "DONE").write_text("ok", encoding="utf-8")
        return record


def one(directory: Path, pattern: str) -> Path:
    hits = sorted(Path(directory).glob(pattern))
    if len(hits) != 1:
        raise RuntimeError(f"expected one {pattern} in {directory}, found {len(hits)}")
    return hits[0]


# --------------------------------------------------------------------------- driver
def job_dir(root: Path, name: str, value, crop: str) -> Path:
    return root / f"{name}_{fmt_value(value)}" / crop


def run_jobs(jobs, root: Path, workers: int, ram_gb: float, keep_tifs: bool) -> None:
    budget = RamBudget(ram_gb)
    mon = PeakMonitor()
    mon.start()
    t0 = time.time()

    def work(job):
        crop, param, direction, value = job
        out = job_dir(root, "identity" if param is None else param["name"], "frozen" if param is None else value, crop)
        if FROZEN in out.parents:
            raise RuntimeError("refusing to write under the frozen runs")
        if (out / "DONE").exists():
            return "skip", job, None
        if out.exists():
            shutil.rmtree(out)
        try:
            CropRun(crop, param, direction, value, out, budget, keep_tifs).execute()
            return "ok", job, None
        except Exception as exc:  # noqa: BLE001 — recorded, the sweep continues
            out.mkdir(parents=True, exist_ok=True)
            (out / "FAILED").write_text(traceback.format_exc(), encoding="utf-8")
            return "failed", job, repr(exc)

    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(work, j) for j in jobs]
        for fut in as_completed(futures):
            status, (crop, param, direction, value), err = fut.result()
            done += 1
            name = "identity" if param is None else param["name"]
            print(f"[sens] {done}/{len(jobs)} {status:6s} {name}={value} {crop} "
                  f"({time.time() - t0:.0f}s, sys peak {(mon.peak - mon.base) / 2**30:.1f} GB above start)"
                  + (f" {err}" if err else ""), flush=True)
    mon.stop.set()
    (root / "_ram_peaks.log").open("a", encoding="utf-8").write(
        f"{dt.datetime.now().isoformat(timespec='seconds')} jobs={len(jobs)} workers={workers} budget={ram_gb}GB "
        f"system_used_peak_above_start={(mon.peak - mon.base) / 2**30:.2f}GB absolute_peak={mon.peak / 2**30:.2f}GB\n")


def collect(root: Path, params: list[dict], out_csv: Path) -> pd.DataFrame:
    """crop-macro mean / min F1@0.50 per perturbation and route; deltas vs the frozen Table II values."""
    base = []
    for crop in CROPS:
        f = FROZEN / crop
        # frozen reference: evaluate_v2 on the frozen segmentations (== e09b seed_ignore, Phase 1d)
        base.append(crop)
    ref = pd.read_csv(REPO_ROOT / "results" / "evaluate_v2" / "per_crop_all_modes.csv")
    ref = ref[(ref["mode"] == "markers_dominant") & (ref["size_filter"].astype(str) == "none")]
    ref_stats = {r: (g["f1_iou050"].mean(), g["f1_iou050"].min(), g["f1_iou075"].mean())
                 for r, g in ref.groupby("route")}
    rows = []
    for p in [None, *params]:
        pairs = [("identity", "frozen")] if p is None else perturbed_values(p)
        for direction, value in pairs:
            name = "identity" if p is None else p["name"]
            per, failed, missing = [], [], []
            for crop in CROPS:
                d = job_dir(root, name, value, crop)
                if (d / "DONE").exists():
                    per.append(pd.read_csv(d / "metrics.csv"))
                elif (d / "FAILED").exists():
                    failed.append(crop)
                else:
                    missing.append(crop)
            if not per and not failed:
                continue
            frame = pd.concat(per, ignore_index=True) if per else pd.DataFrame(columns=["route", "f1_iou050"])
            for route in ("BA-TOS_sel", "BA-TOS_ToS"):
                g = frame[frame["route"] == route]
                complete = len(g) == len(CROPS)
                mean = g["f1_iou050"].mean() if complete else np.nan
                mn = g["f1_iou050"].min() if complete else np.nan
                rm, rmin, r75 = ref_stats[route]
                row = {"group": "identity" if p is None else p["group"], "param": name,
                       "kind": "identity" if p is None else p["kind"], "stage": "" if p is None else p.get("stage", ""),
                       "script": "" if p is None else p.get("script", ""),
                       "frozen_value": "" if p is None else p["frozen"], "direction": direction,
                       "perturbed_value": value, "route": route, "n_crops_done": len(g),
                       "failed_crops": ";".join(failed), "missing_crops": ";".join(missing),
                       "mean_f1_050": mean, "min_f1_050": mn,
                       "mean_f1_075": g["f1_iou075"].mean() if complete else np.nan,
                       "delta_mean_f1_050": mean - rm, "delta_min_f1_050": mn - rmin,
                       "n_gate_triggered": int(g["gate_triggered"].sum()) if complete and "gate_triggered" in g else np.nan}
                rows.append(row)
    table = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_csv, index=False)
    return table


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Phase 5a global sensitivity")
    ap.add_argument("--config", type=Path, default=REPO_ROOT / "configs" / "sensitivity_global.yaml")
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--groups", nargs="*", help="groups to run (default: all)")
    ap.add_argument("--params", nargs="*", help="parameter names to run (default: all in the groups)")
    ap.add_argument("--identity", action="store_true", help="also run the unperturbed pipeline from scratch")
    ap.add_argument("--values", nargs="*",
                    help="only these perturbed values (matched on the value as written, e.g. 2.0, or its folder tag, e.g. 2)")
    ap.add_argument("--crops", nargs="*", default=CROPS)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--ram-gb", type=float, default=22.0)
    ap.add_argument("--keep-tifs", action="store_true")
    ap.add_argument("--collect-only", action="store_true")
    ap.add_argument("--out-csv", type=Path, default=REPO_ROOT / "results" / "sensitivity" / "global_sensitivity.csv")
    args = ap.parse_args(argv)
    params = load_params(args.config)
    selected = [p for p in params if (not args.groups or p["group"] in args.groups)
                and (not args.params or p["name"] in args.params)]
    if not args.collect_only:
        jobs = []
        if args.identity:
            jobs += [(c, None, "identity", "frozen") for c in args.crops]
        for p in selected:
            for direction, value in perturbed_values(p):
                if args.values and str(value) not in args.values and fmt_value(value) not in args.values:
                    continue
                jobs += [(c, p, direction, value) for c in args.crops]
        # heavy-first ordering keeps the RAM budget busy
        args.root.mkdir(parents=True, exist_ok=True)
        (args.root / "_git.json").write_text(json.dumps(git_state(), indent=1), encoding="utf-8")
        print(f"[sens] {len(jobs)} crop-jobs, workers={args.workers}, RAM budget {args.ram_gb} GB", flush=True)
        run_jobs(jobs, args.root, args.workers, args.ram_gb, args.keep_tifs)
    table = collect(args.root, params, args.out_csv)
    print(table.to_string(index=False, max_colwidth=30))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
