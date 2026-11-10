#!/usr/bin/env python3
"""Phase 5i (protocol clarification): synthetic scale crops with a canonical 200^3 field of view.

For each dev crop and synthetic factor F in {0.7, 1.4, 2.0} (grains scaled by F, as in make_synthetic_scale.py):
  probe     = the existing 200^3 synthetic crop at the same centre (crops_synth_scale.json);
              r_probe = its r_granulo (batos_v3.scale), recomputed here and checked against the 5e record
  native S  = round(200 * r_probe / 10)          (so that S * 10 / r_probe = 200 canonical voxels)
  window    = the full EFRGP01_00/01 volumes resampled by F about the dev crop centre, central S^3:
              output voxel o <- source c + (o - (S - 1) / 2) / F   (grey: linear, rounded uint16; labels: nearest)
              (identical transform to make_synthetic_scale.py; S = 200 reproduces the 5e crop exactly: checked)
  valid GT  = >= round(1000 * (r_probe / 10)^3) native voxels (= 1000 at canonical scale), not touching a face.
GT is not used to define any window. Outputs: C:\\BA_TOS_DATA\\crops_synth_canon\\x<F>\\<crop_id>\\{gray,labels}.tif,
data/metadata/crops_synth_canon.json, data/metadata/crops_manifest_synth_canon.csv.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.evaluate_v2 import reference_protocol  # noqa: E402
from batos_v3.scale import estimate  # noqa: E402

DEV = REPO / "data" / "metadata" / "crops_sibgrapi.json"
SYN = REPO / "data" / "metadata" / "crops_synth_scale.json"
V3INFO = Path(r"C:\Users\USER\batos-3d-grain-segmentation\results\v3_scale_test\crop_info.csv")
OUT_JSON = REPO / "data" / "metadata" / "crops_synth_canon.json"
OUT_CSV = REPO / "data" / "metadata" / "crops_manifest_synth_canon.csv"
ROOT = Path(r"C:\BA_TOS_DATA\crops_synth_canon")
R_REF, CANON = 10.0, 200


def load_make_crops():
    spec = importlib.util.spec_from_file_location("make_crops_ro", REPO / "scripts" / "make_crops.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tag(f: float) -> str:
    return f"x{f:g}".replace(".", "p")


def resample(vol, centre, f: float, size: int, order: int, dtype) -> tuple[np.ndarray, list]:
    half = (size - 1) / 2.0
    lo, hi = [], []
    for c, n in zip(centre, vol.shape):
        lo.append(int(math.floor(c - half / f)) - 1)
        hi.append(int(math.ceil(c + half / f)) + 2)
        if lo[-1] < 0 or hi[-1] > n:
            raise RuntimeError(f"source window {lo[-1]}..{hi[-1]} outside volume extent {n}")
    block = np.asarray(vol[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]])
    offset = [c - l - half / f for c, l in zip(centre, lo)]
    work = block.astype(np.float32) if order else block
    out = ndi.affine_transform(work, matrix=[1.0 / f] * 3, offset=offset, output_shape=(size,) * 3, order=order,
                               mode="nearest", prefilter=False)
    if order:
        out = np.clip(np.rint(out), np.iinfo(dtype).min, np.iinfo(dtype).max)
    return out.astype(dtype), [lo, hi]


def sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def main() -> int:
    import pandas as pd
    mc = load_make_crops()
    dev = json.loads(DEV.read_text(encoding="utf-8"))
    syn = {c["crop_id"]: c for c in json.loads(SYN.read_text(encoding="utf-8"))["crops"]}
    rec5e = pd.read_csv(V3INFO).set_index("crop_id").v3_r_scale.to_dict()
    crops_out, rows = [], []
    for scan, sv in dev["source_volumes"].items():
        if not isinstance(sv, dict) or "gray_path" not in sv:
            continue
        gray_v = mc.open_memmap(Path(sv["gray_path"]))
        gt_v = mc.open_memmap(Path(sv["gt_full_path"]))
        for c in (c for c in dev["crops"] if c["scan"] == scan):
            centre = [c["z"] + 99.5, c["y"] + 99.5, c["x"] + 99.5]
            for F in (0.7, 1.4, 2.0):
                pid = f"{c['crop_id']}_{tag(F)}"
                probe = tifffile.imread(syn[pid]["gray_crop_path"])
                p200, _ = resample(gray_v, centre, F, CANON, 1, gray_v.dtype)
                if sha(p200) != sha(probe):
                    raise RuntimeError(f"S = 200 resample differs from the 5e synthetic crop {pid}")
                r = float(estimate(probe)[0]["r_granulo"])
                if r != float(rec5e[pid]):
                    raise RuntimeError(f"{pid}: r_probe {r} != recorded 5e r_granulo {rec5e[pid]}")
                S = int(round(CANON * r / R_REF))
                cid = f"{pid}_canon"
                d = ROOT / tag(F) / cid
                d.mkdir(parents=True, exist_ok=True)
                g, box = resample(gray_v, centre, F, S, 1, gray_v.dtype)
                lab, _ = resample(gt_v, centre, F, S, 0, gt_v.dtype)
                tifffile.imwrite(d / "gray.tif", g)
                tifffile.imwrite(d / "labels.tif", lab)
                del g
                mgv = int(round(1000 * (r / R_REF) ** 3))
                n_valid = int(len(reference_protocol(lab, mgv, True)["valid_ids"]))
                crops_out.append({"crop_id": cid, "dev_crop_id": c["crop_id"], "probe_crop_id": pid, "scan": scan,
                                  "specimen": "EFR", "factor": F, "r_probe": r, "native_side": S,
                                  "canonical_side": CANON, "centre_zyx": centre, "source_box_zyx": box,
                                  "gray_crop_path": str(d / "gray.tif"), "gt_path": str(d / "labels.tif"),
                                  "min_gt_voxels": mgv})
                rows.append({"crop_id": cid, "factor": F, "r_probe": r, "native_side": S,
                             "gray_sha256": sha(tifffile.imread(d / "gray.tif")), "labels_sha256": sha(lab),
                             "min_gt_voxels": mgv, "valid_gt": n_valid})
                print(f"[canon] {cid}: r_probe {r:g} -> native {S}^3, valid GT {n_valid} (min {mgv} vox)", flush=True)
                del lab
    meta = {"schema_version": 1, "created": "2026-11-10", "description": __doc__, "factors": [0.7, 1.4, 2.0],
            "crop_convention": "zyx", "crops": crops_out}
    OUT_JSON.write_text(json.dumps(meta, indent=1), encoding="utf-8")
    with OUT_CSV.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"[done] {len(crops_out)} crops -> {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
