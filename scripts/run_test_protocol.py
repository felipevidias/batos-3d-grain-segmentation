#!/usr/bin/env python3
"""Phase 1f: frozen SIBGRAPI pipeline (with gate) and batos-paper-v2-fixflat on the 52 test crops, then aggregation.

No parameter is tuned; both methods are run exactly as frozen (configs/batos_paper_v2_fixflat.yaml).
Per crop, sequentially (own peak = one stage, < 8 GB):
  0. e13 scale (verbatim) per specimen. e13 requires exactly 12 manifest rows; specimens with 8 test crops are padded
     with 4 duplicate rows (crop_id + "__pad"), whose outputs are discarded. e13 is per-crop independent (each row is
     estimated from its own crop only), so the padding does not change any estimate.
  1. frozen: python -I -m batos_v2.frozen --crops-json data/metadata/crops_test.json --openrec-impl auto
     (openrec radius > 15 -> batos_v2.openrec_lowmem, bit-identical to the frozen script; recorded per crop);
  2. fixflat: batos_v2.flatten_fix_run.FlattenFixRun variant A, upstream (scale, openrec) read from the step-1 run;
  3. evaluate_v2 on both methods x routes (BA-TOS_ToS, BA-TOS_sel) x size filter (none, 500: labels < 500 voxels
     zeroed, no relabel, as filter_predicted_labels_by_size.py) x mode (markers_dominant = e09b seed-ignore;
     no_markers_fraction050 = own-voxel support, overlap-ignore 0.50), IoU 0.25 / 0.50 / 0.75;
  4. scale diagnostics: e13 r_est, granulometric r (batos_v2.scale_granulo: r_ps_peak, r_ps_median), GT median
     inscribed radius over valid GT grains (GT for diagnosis only);
  5. prune intermediate TIFFs (disk budget): keep the ToS / sel segmentations and their normalized markers.
--aggregate writes results/test_protocol/*.csv, REPORT tables and results/all_runs.csv rows.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2 import sensitivity  # noqa: E402
from batos_v2.evaluate_v2 import evaluate  # noqa: E402
from batos_v2.extinction_run import upsert  # noqa: E402
from batos_v2.flatten_fix_run import FlattenFixRun  # noqa: E402
from batos_v2.frozen import SCRIPTS  # noqa: E402
from batos_v2.run_crop import git_state  # noqa: E402
from batos_v2.scale_granulo import estimate as granulo_estimate  # noqa: E402

PY = sys.executable
META = REPO / "data" / "metadata" / "crops_test.json"
ROOT = Path(r"C:\BA_TOS_DATA\runs\test_protocol")
FROZEN_ROOT = ROOT / "frozen"
FIX_ROOT = ROOT / "fixflat"
OUT = REPO / "results" / "test_protocol"
ROUTES = {"BA-TOS_ToS": ("e11/watershed", "tos_e11_min"), "BA-TOS_sel": ("e23c/watershed", "batos_selective")}
MODES = {"markers_dominant": dict(rule="dominant"),
         "no_markers_fraction050": dict(rule="fraction", ignore_overlap_threshold=0.50)}
KEEP = ("segmentations/tos_e11_min_distance_c1.tif", "segmentations/batos_selective_distance_c1.tif",
        "normalized_markers/tos_e11_min.tif", "normalized_markers/batos_selective.tif")


def log(msg: str) -> None:
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def scale_stage(meta: dict, spec: str) -> None:
    """e13 (verbatim) for one specimen's test crops; outputs merged into <FROZEN_ROOT>/_scale/auto_scale_params."""
    params = FROZEN_ROOT / "_scale" / "auto_scale_params"
    params.mkdir(parents=True, exist_ok=True)
    d = FROZEN_ROOT / "_scale" / spec
    crops = [c for c in meta["crops"] if c["specimen"] == spec]
    if not (d / "out" / "e13_auto_scale_params_per_crop.csv").exists():
        rows = [{"scan": c["scan"], "crop_id": c["crop_id"], "gray": c["gray_crop_path"]} for c in crops]
        pads = [{**r, "crop_id": r["crop_id"] + "__pad"} for r in rows[:12 - len(rows)]]
        d.mkdir(parents=True, exist_ok=True)
        man = d / "scale_manifest.csv"
        with man.open("w", encoding="utf-8", newline="") as h:
            w = csv.DictWriter(h, fieldnames=["scan", "crop_id", "gray"])
            w.writeheader()
            w.writerows(rows + pads)
        t0 = time.time()
        subprocess.run([PY, "-I", str(SCRIPTS["e13_estimate_scale"]), "--manifest", str(man), "--out-dir",
                        str(d / "out")], check=True, stdout=subprocess.DEVNULL)
        per = pd.read_csv(d / "out" / "e13_auto_scale_params_per_crop.csv")
        for pad in pads:  # padding must reproduce the original row exactly (per-crop independence check)
            a = per[per.crop_id == pad["crop_id"]].drop(columns="crop_id").iloc[0]
            b = per[per.crop_id == pad["crop_id"][:-5]].drop(columns="crop_id").iloc[0]
            if not a.astype(str).equals(b.astype(str)):
                raise RuntimeError(f"e13 padding row differs from original: {pad['crop_id']}")
        log(f"scale {spec}: {len(crops)} crops (+{len(pads)} pad rows) in {time.time() - t0:.0f}s; r_est "
            f"{per[~per.crop_id.str.endswith('__pad')].estimated_grain_radius.round(2).tolist()}")
    for c in crops:
        shutil.copy2(d / "out" / f"{c['crop_id']}_auto_params.json", params / f"{c['crop_id']}_auto_params.json")
    merged = [pd.read_csv(f) for f in sorted((FROZEN_ROOT / "_scale").glob("*/out/e13_auto_scale_params_per_crop.csv"))]
    merged = pd.concat(merged)
    merged[~merged.crop_id.str.endswith("__pad")].to_csv(params / "e13_auto_scale_params_per_crop.csv", index=False)


def run_frozen(crop: str) -> None:
    out = FROZEN_ROOT / crop
    if (out / "frozen_run.json").exists():
        return
    if out.exists():
        shutil.rmtree(out)  # partial test run (our own output, never the dev frozen runs)
    t0 = time.time()
    with open(ROOT / "logs" / f"{crop}_frozen.log", "w", encoding="utf-8") as h:
        rc = subprocess.run([PY, "-I", "-m", "batos_v2.frozen", "--crop-id", crop, "--crops-json", str(META),
                             "--runs-root", str(FROZEN_ROOT), "--openrec-impl", "auto"],
                            stdout=h, stderr=subprocess.STDOUT, cwd=REPO / "src").returncode
    if rc != 0:
        raise RuntimeError(f"frozen failed on {crop} (rc {rc}); see {ROOT / 'logs' / f'{crop}_frozen.log'}")
    log(f"{crop} frozen done in {time.time() - t0:.0f}s")


def run_fixflat(crop: str, git: dict) -> None:
    out = FIX_ROOT / crop
    if (out / "DONE").exists():
        return
    if out.exists():
        shutil.rmtree(out)
    sensitivity.FROZEN = FROZEN_ROOT  # CropRun.frozen_paths reads <FROZEN>/<crop>/frozen_run.json
    param = {"name": "batos_paper_v2_fixflat", "group": "test_protocol", "kind": "flag", "stage": "flatten",
             "flag": "--variant", "frozen": "frozen", "type": "str"}
    t0 = time.time()
    run = FlattenFixRun(crop, param, "fix", "A", out, sensitivity.RamBudget(8.0), keep_tifs=True)
    run.note = "flatten stage = batos_v2.flatten_fixed --variant A; everything else frozen (test protocol)"
    run.execute()
    rec = json.loads((out / "run.json").read_text(encoding="utf-8"))
    info = json.loads((out / "pipe" / "01_preprocess" / "flatten_fixed_info.json").read_text(encoding="utf-8"))
    rec |= {"git": git, "bg_ref": info["bg_ref"]}
    (out / "run.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    log(f"{crop} fixflat done in {time.time() - t0:.0f}s (bg_ref {info['bg_ref']:.1f})")


def size_filter(pred: np.ndarray, min_size: int = 500) -> np.ndarray:
    ids, counts = np.unique(pred, return_counts=True)
    keep = ids[(ids > 0) & (counts >= min_size)]
    return np.where(np.isin(pred, keep), pred, 0).astype(pred.dtype)


def inscribed_radii(gt: np.ndarray) -> np.ndarray:
    """Max EDT inside each valid GT grain (>= 1000 voxels, not touching a face). GT diagnostic only."""
    faces = np.unique(np.concatenate([gt[0].ravel(), gt[-1].ravel(), gt[:, 0].ravel(), gt[:, -1].ravel(),
                                      gt[:, :, 0].ravel(), gt[:, :, -1].ravel()]))
    out = []
    for lab, sl in enumerate(ndi.find_objects(gt), start=0):
        if sl is None:
            continue
        lab_id = lab + 1
        if lab_id in faces:
            continue
        m = gt[sl] == lab_id
        if m.sum() < 1000:
            continue
        out.append(ndi.distance_transform_edt(np.pad(m, 1)).max())
    return np.asarray(out)


def evaluate_crop(crop: str, c: dict) -> None:
    done = ROOT / "eval" / f"{crop}.json"
    if done.exists():
        return
    t0 = time.time()
    gt = tifffile.imread(c["gt_path"])
    rows = []
    for method, base in (("frozen", FROZEN_ROOT / crop), ("fixflat", FIX_ROOT / crop)):
        for route, (ws, label) in ROUTES.items():
            pred = tifffile.imread(base / ws / "segmentations" / f"{label}_distance_c1.tif")
            markers = tifffile.imread(base / ws / "normalized_markers" / f"{label}.tif")
            for flt in ("none", "500"):
                p = pred if flt == "none" else size_filter(pred)
                for mode, kw in MODES.items():
                    m, _, _ = evaluate(p, gt, markers if mode.startswith("markers") else None, morph_iou=None, **kw)
                    rows.append({"crop_id": crop, "specimen": c["specimen"], "scan": c["scan"],
                                 "test_class": c["test_class"], "method": method, "route": route,
                                 "size_filter": flt, "mode": mode, **m})
    frec = json.loads((FROZEN_ROOT / crop / "frozen_run.json").read_text(encoding="utf-8"))
    xrec = json.loads((FIX_ROOT / crop / "run.json").read_text(encoding="utf-8"))
    gate = {}
    for method, rec in (("frozen", frec), ("fixflat", xrec)):
        g = rec["route"].get("gate", rec["route"]) if method == "frozen" else rec["gate"]
        gate[method] = {"gate_triggered": bool(g["gate_triggered"]), "mv": g["region_voxels_median"],
                        "qv": g["region_voxels_q75"]}
    n_markers = {}
    for method, base in (("frozen", FROZEN_ROOT / crop), ("fixflat", FIX_ROOT / crop)):
        mk = pd.read_csv(base / "e11" / "markers_tos_min" / "markers.csv") \
            if (base / "e11" / "markers_tos_min" / "markers.csv").exists() else None
        n_markers[method] = None if mk is None else int(len(mk))
    gran = granulo_estimate(tifffile.imread(c["gray_crop_path"]))
    rin = inscribed_radii(gt)
    info = {"crop_id": crop, "specimen": c["specimen"], "scan": c["scan"], "test_class": c["test_class"],
            "specimen_fraction_sub4": c["specimen_fraction_sub4"],
            "r_est_e13": frec["auto_params"]["estimated_grain_radius"],
            "openrec_radius": frec["openrec_radius"], "openrec_implementation": frec.get("openrec_implementation"),
            "r_ps_peak": gran.get("r_ps_peak"), "r_ps_median": gran.get("r_ps_median"),
            "gt_inscribed_r_median": float(np.median(rin)) if rin.size else None, "gt_valid_n": int(rin.size),
            "bg_ref_fixflat": xrec.get("bg_ref"),
            **{f"{m}_{k}": v for m, g in gate.items() for k, v in g.items()},
            "frozen_tos_markers": n_markers["frozen"], "fixflat_tos_markers": n_markers["fixflat"],
            "frozen_seconds": frec["seconds_total"], "frozen_peak_mb": frec.get("peak_memory_mb"),
            "fixflat_seconds_flatten_onward": round(sum(s["seconds"] for s in xrec["steps"]), 1),
            "fixflat_peak_mb": max(s.get("peak_memory_mb", 0) for s in xrec["steps"]),
            "eval_seconds": None}
    upstream = sum(s["seconds"] for s in frec["steps"] if s["step"] in ("scale", "openrec_v2"))
    info["fixflat_seconds_end_to_end"] = round(info["fixflat_seconds_flatten_onward"] + upstream, 1)
    info["eval_seconds"] = round(time.time() - t0, 1)
    (ROOT / "eval").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(ROOT / "eval" / f"{crop}.csv", index=False)
    done.write_text(json.dumps(info, indent=1, default=str), encoding="utf-8")
    f = pd.DataFrame(rows)
    sel = f[(f["mode"] == "markers_dominant") & (f.size_filter == "none")]
    log(f"{crop} eval: " + ", ".join(f"{r.method}/{r.route[7:]} {r.f1_iou050:.3f}" for r in sel.itertuples())
        + f" | gate frozen {gate['frozen']} fixflat {gate['fixflat']}")


def prune(crop: str) -> None:
    for base in (FROZEN_ROOT / crop, FIX_ROOT / crop):
        for p in base.rglob("*.tif"):
            rel = p.relative_to(base).as_posix()
            if not any(rel.endswith(k) for k in KEEP):
                p.unlink()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--crops", nargs="*")
    ap.add_argument("--no-prune", action="store_true")
    args = ap.parse_args(argv)
    meta = json.loads(META.read_text(encoding="utf-8"))
    crops = [c for c in meta["crops"] if not args.crops or c["crop_id"] in args.crops]
    for d in (FROZEN_ROOT, FIX_ROOT, ROOT / "logs", ROOT / "eval"):
        d.mkdir(parents=True, exist_ok=True)
    git = git_state()
    for i, c in enumerate(crops, 1):
        crop = c["crop_id"]
        if not Path(c["gt_path"]).exists():
            raise FileNotFoundError(f"crop not cut yet: {c['gt_path']}")
        log(f"--- {i}/{len(crops)} {crop} ({c['test_class']})")
        scale_stage(meta, c["specimen"])
        run_frozen(crop)
        run_fixflat(crop, git)
        evaluate_crop(crop, c)
        if not args.no_prune:
            prune(crop)
    log("END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
