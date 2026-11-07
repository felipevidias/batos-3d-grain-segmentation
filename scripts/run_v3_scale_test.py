#!/usr/bin/env python3
"""Phase 5e step 3: synthetic scale test, fixflat vs v3, on the 12 EFR dev crops at f = 1.0 and their rescaled copies
at f = 0.7 / 1.4 / 2.0 (data/metadata/crops_synth_scale.json). Dev + synthetic only; no test crop is read.

Per factor f != 1: e13 (verbatim, the 12 crops of that factor = exactly 12 manifest rows) -> fixflat preset; v3 preset.
f = 1.0: v3 runs here; fixflat = the existing batos-paper-v2-fixflat segmentations (C:\\BA_TOS_DATA\\runs\\
batos_paper_v2_fixflat\\A), which the batos_v3.pipeline fixflat preset reproduces bit-for-bit (checked on crop04).
Evaluation: batos_v2.evaluate_v2, routes ToS / sel, modes markers_dominant (seed-ignore, primary) and
no_markers_fraction050, size filter off and the method's own (fixflat 500; v3 500 s^3), IoU .25/.50/.75, valid GT
>= round(1000 f^3) voxels (same physical grains at every factor). GT inscribed radius = diagnostic only.
Outputs: C:\\BA_TOS_DATA\\runs\\v3_scale_test\\{e13,runs,eval}; log C:\\BA_TOS_DATA\\runs\\v3_scale_test.log.
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
from batos_v2.evaluate_v2 import evaluate  # noqa: E402
from batos_v2.frozen import SCRIPTS  # noqa: E402
from batos_v3.pipeline import StageError, outputs, run_crop  # noqa: E402

PY = sys.executable
ROOT = Path(r"C:\BA_TOS_DATA\runs\v3_scale_test")
DEV_FIXFLAT = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_fixflat\A")
DEV_E13 = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen\_scale\auto_scale_params")
MODES = {"markers_dominant": dict(rule="dominant"),
         "no_markers_fraction050": dict(rule="fraction", ignore_overlap_threshold=0.50)}
KEEP = ("segmentations/tos_e11_min_distance_c1.tif", "segmentations/batos_selective_distance_c1.tif",
        "normalized_markers/tos_e11_min.tif", "normalized_markers/batos_selective.tif")


def log(msg: str) -> None:
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def crop_sets() -> dict[float, list[dict]]:
    dev = json.loads((REPO / "data/metadata/crops_sibgrapi.json").read_text(encoding="utf-8"))
    syn = json.loads((REPO / "data/metadata/crops_synth_scale.json").read_text(encoding="utf-8"))
    sets = {1.0: [{"crop_id": c["crop_id"], "dev_crop_id": c["crop_id"], "scan": c["scan"], "factor": 1.0,
                   "gray_crop_path": c["gray_crop_path"], "gt_path": c["gt_path"], "min_gt_voxels": 1000}
                  for c in dev["crops"]]}
    for c in syn["crops"]:
        sets.setdefault(float(c["factor"]), []).append(c)
    return sets


def e13_for(f: float, crops: list[dict]) -> Path:
    if f == 1.0:
        return DEV_E13
    d = ROOT / "e13" / f"x{f:g}"
    out = d / "out"
    if not (out / "e13_auto_scale_params_per_crop.csv").exists():
        d.mkdir(parents=True, exist_ok=True)
        man = d / "scale_manifest.csv"
        with man.open("w", encoding="utf-8", newline="") as h:
            w = csv.DictWriter(h, fieldnames=["scan", "crop_id", "gray"])
            w.writeheader()
            w.writerows({"scan": c["scan"], "crop_id": c["crop_id"], "gray": c["gray_crop_path"]} for c in crops)
        subprocess.run([PY, "-I", str(SCRIPTS["e13_estimate_scale"]), "--manifest", str(man), "--out-dir", str(out)],
                       check=True, stdout=subprocess.DEVNULL)
        per = pd.read_csv(out / "e13_auto_scale_params_per_crop.csv")
        log(f"e13 x{f:g}: r_est {per.estimated_grain_radius.round(2).tolist()}")
    return out


def run_dir(method: str, f: float, crop: str) -> Path:
    if method == "fixflat" and f == 1.0:
        return DEV_FIXFLAT / crop
    return ROOT / "runs" / method / f"x{f:g}" / crop


def ensure_run(method: str, f: float, c: dict, e13_dir: Path) -> None:
    out = run_dir(method, f, c["crop_id"])
    if method == "fixflat" and f == 1.0:
        return
    if (out / "DONE").exists() or (out / "EMPTY").exists() or (out / "FAILED.json").exists():
        return
    if out.exists():
        shutil.rmtree(out)
    t0 = time.time()
    try:
        rec = run_crop(c["crop_id"], Path(c["gray_crop_path"]), out, method,
                       e13_dir / f"{c['crop_id']}_auto_params.json" if method == "fixflat" else None)
        log(f"{c['crop_id']} {method} done in {time.time() - t0:.0f}s | r {rec.get('r_granulo', rec.get('r_est_e13'))} "
            f"openrec {rec['openrec_radius']} thr03 {rec.get('stage03_dark_threshold')} gate {rec.get('gate')} "
            f"markers {rec.get('tos_markers')} {'EMPTY' if rec.get('empty_result') else ''}")
    except StageError as e:  # recorded method failure, evaluated as an empty segmentation
        out.mkdir(parents=True, exist_ok=True)
        (out / "FAILED.json").write_text(json.dumps({"error": str(e), "seconds": round(time.time() - t0, 1)}),
                                         encoding="utf-8")
        log(f"{c['crop_id']} {method} FAILED: {e}")


def inscribed_median(gt: np.ndarray, min_vox: int) -> float | None:
    faces = set(np.unique(np.concatenate([gt[0].ravel(), gt[-1].ravel(), gt[:, 0].ravel(), gt[:, -1].ravel(),
                                          gt[:, :, 0].ravel(), gt[:, :, -1].ravel()])).tolist())
    vals = []
    for i, sl in enumerate(ndi.find_objects(gt)):
        if sl is None or (i + 1) in faces:
            continue
        m = gt[sl] == i + 1
        if m.sum() >= min_vox:
            vals.append(ndi.distance_transform_edt(np.pad(m, 1)).max())
    return float(np.median(vals)) if vals else None


def run_record(method: str, f: float, crop: str) -> dict:
    out = run_dir(method, f, crop)
    if method == "fixflat" and f == 1.0:
        rec = json.loads((out / "run.json").read_text(encoding="utf-8"))
        g = rec["gate"]
        return {"gate": {"gate_triggered": g["gate_triggered"], "region_voxels_median": g["region_voxels_median"],
                         "region_voxels_q75": g["region_voxels_q75"]}, "openrec_radius": rec["openrec_radius"],
                "seconds_total": None, "peak_memory_mb": max(s.get("peak_memory_mb", 0) for s in rec["steps"]),
                "size_filter": 500, "note": "existing batos-paper-v2-fixflat run (flatten onward)"}
    for name in ("v3_run.json", "FAILED.json"):
        if (out / name).exists():
            return json.loads((out / name).read_text(encoding="utf-8"))
    raise FileNotFoundError(out)


def size_filter(pred: np.ndarray, min_size: int) -> np.ndarray:
    ids, counts = np.unique(pred, return_counts=True)
    keep = ids[(ids > 0) & (counts >= min_size)]
    return np.where(np.isin(pred, keep), pred, 0).astype(pred.dtype)


def evaluate_crop(f: float, c: dict) -> None:
    done = ROOT / "eval" / f"{c['crop_id']}.json"
    if done.exists():
        return
    gt = tifffile.imread(c["gt_path"])
    rows, info = [], {"crop_id": c["crop_id"], "dev_crop_id": c["dev_crop_id"], "factor": f,
                      "min_gt_voxels": c["min_gt_voxels"],
                      "gt_inscribed_r_median": inscribed_median(gt, c["min_gt_voxels"])}
    for method in ("fixflat", "v3"):
        out = run_dir(method, f, c["crop_id"])
        rec = run_record(method, f, c["crop_id"])
        failed = (out / "FAILED.json").exists()
        res = outputs(out, gt.shape) if not failed else {r: (np.zeros(gt.shape, np.uint32),) * 2
                                                         for r in ("BA-TOS_ToS", "BA-TOS_sel")}
        own = int(rec.get("size_filter", 500))
        for route, (pred, markers) in res.items():
            for flt, ms in (("none", None), ("own", own)):
                p = pred if ms is None else size_filter(pred, ms)
                for mode, kw in MODES.items():
                    m, _, _ = evaluate(p, gt, markers if mode.startswith("markers") else None, morph_iou=None,
                                       min_gt_voxels=c["min_gt_voxels"], **kw)
                    rows.append({"crop_id": c["crop_id"], "dev_crop_id": c["dev_crop_id"], "factor": f,
                                 "method": method, "route": route, "size_filter": flt, "size_filter_voxels": ms,
                                 "mode": mode, **m})
        g = rec.get("gate") or {}
        info |= {f"{method}_r_scale": rec.get("r_granulo", rec.get("r_est_e13")), f"{method}_s": rec.get("s"),
                 f"{method}_openrec_radius": rec.get("openrec_radius"),
                 f"{method}_openrec_impl": rec.get("openrec_implementation"),
                 f"{method}_stage03_thr": rec.get("stage03_dark_threshold"),
                 f"{method}_gate": g.get("gate_triggered"), f"{method}_mv": g.get("region_voxels_median"),
                 f"{method}_qv": g.get("region_voxels_q75"), f"{method}_tos_markers": rec.get("tos_markers"),
                 f"{method}_empty": rec.get("empty_result"), f"{method}_failed": rec.get("error"),
                 f"{method}_seconds": rec.get("seconds_total"), f"{method}_peak_mb": rec.get("peak_memory_mb"),
                 f"{method}_size_filter": own}
    if info["fixflat_r_scale"] is None and f == 1.0:
        info["fixflat_r_scale"] = json.loads((DEV_E13 / f"{c['crop_id']}_auto_params.json").read_text(
            encoding="utf-8"))["estimated_grain_radius"]
    (ROOT / "eval").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(ROOT / "eval" / f"{c['crop_id']}.csv", index=False)
    done.write_text(json.dumps(info, indent=1, default=str), encoding="utf-8")
    sel = pd.DataFrame(rows)
    sel = sel[(sel["mode"] == "markers_dominant") & (sel.size_filter == "none")]
    log(f"{c['crop_id']} eval: " + ", ".join(f"{r.method}/{r.route[7:]} {r.f1_iou050:.3f}" for r in sel.itertuples()))


def prune(f: float, crop: str) -> None:
    for method in ("fixflat", "v3"):
        if method == "fixflat" and f == 1.0:
            continue
        for p in run_dir(method, f, crop).rglob("*.tif"):
            if not any(p.as_posix().endswith(k) for k in KEEP):
                p.unlink()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--factors", nargs="*", type=float, default=[1.0, 1.4, 2.0, 0.7])
    ap.add_argument("--crops", nargs="*", help="dev crop ids (default all 12)")
    ap.add_argument("--no-prune", action="store_true")
    a = ap.parse_args(argv)
    sets = crop_sets()
    for f in a.factors:
        crops = [c for c in sets[f] if not a.crops or c["dev_crop_id"] in a.crops]
        e13_dir = e13_for(f, sets[f])
        for i, c in enumerate(crops, 1):
            log(f"--- x{f:g} {i}/{len(crops)} {c['crop_id']}")
            for method in ("fixflat", "v3"):
                ensure_run(method, f, c, e13_dir)
            evaluate_crop(f, c)
            if not a.no_prune:
                prune(f, c["crop_id"])
    log("END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
