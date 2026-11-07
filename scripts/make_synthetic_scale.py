#!/usr/bin/env python3
"""Phase 5e step 3: synthetic scale test crops from the 12 EFR dev crops (never from test crops).

For each dev crop (origin x/y/z, size 200, zyx) and factor f in {0.7, 1.4, 2.0}: the full EFRGP01_00/01 volumes
(greyscale and GT labels) are resampled by f about the dev crop centre c = origin + 99.5 and the central 200^3 is kept:
  output voxel o (per axis) <- source coordinate c + (o - 99.5) / f
  greyscale: linear interpolation (scipy affine_transform order 1), rounded to uint16
  labels   : nearest neighbour (order 0), labels unchanged
f = 1 reproduces the dev crop exactly (integer shift), so the 1.0 point uses the original dev crops.
For f < 1 the source window (286^3 at 0.7) extends beyond the dev crop; all windows lie inside the volume.
Evaluation threshold for valid GT scales with the voxel volume: min_gt_voxels = round(1000 f^3) (same physical grains).
Outputs: C:\\BA_TOS_DATA\\crops_synth\\x<f>\\<crop_id>\\{gray,labels}.tif, data/metadata/crops_synth_scale.json,
data/metadata/crops_manifest_synth.csv.
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

DEV = REPO / "data" / "metadata" / "crops_sibgrapi.json"
OUT_JSON = REPO / "data" / "metadata" / "crops_synth_scale.json"
OUT_CSV = REPO / "data" / "metadata" / "crops_manifest_synth.csv"
ROOT = Path(r"C:\BA_TOS_DATA\crops_synth")
FACTORS = (0.7, 1.4, 2.0)
SIZE = 200


def load_make_crops():
    spec = importlib.util.spec_from_file_location("make_crops_ro", REPO / "scripts" / "make_crops.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tag(f: float) -> str:
    return f"x{f:g}".replace(".", "p")


def resample(vol, centre, f: float, order: int, dtype) -> tuple[np.ndarray, list]:
    half = (SIZE - 1) / 2.0  # 99.5
    lo, hi = [], []
    for c, n in zip(centre, vol.shape):
        a, b = c - half / f, c + half / f
        lo.append(int(math.floor(a)) - 1)
        hi.append(int(math.ceil(b)) + 2)
        if lo[-1] < 0 or hi[-1] > n:
            raise RuntimeError(f"source window {lo[-1]}..{hi[-1]} outside volume extent {n}")
    block = np.asarray(vol[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]])
    offset = [c - l - half / f for c, l in zip(centre, lo)]
    work = block.astype(np.float32) if order else block
    out = ndi.affine_transform(work, matrix=[1.0 / f] * 3, offset=offset, output_shape=(SIZE,) * 3, order=order,
                               mode="nearest", prefilter=False)
    if order:
        out = np.clip(np.rint(out), np.iinfo(dtype).min, np.iinfo(dtype).max)
    return out.astype(dtype), [lo, hi]


def sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def main() -> int:
    mc = load_make_crops()
    dev = json.loads(DEV.read_text(encoding="utf-8"))
    crops_out, rows = [], []
    for scan, sv in dev["source_volumes"].items():
        if not isinstance(sv, dict) or "gray_path" not in sv:
            continue
        gray_v = mc.open_memmap(Path(sv["gray_path"]))
        gt_v = mc.open_memmap(Path(sv["gt_full_path"]))
        for c in (c for c in dev["crops"] if c["scan"] == scan):
            centre = [c["z"] + 99.5, c["y"] + 99.5, c["x"] + 99.5]
            # f = 1 identity check against the committed dev crop
            g1, _ = resample(gray_v, centre, 1.0, 1, gray_v.dtype)
            if sha(g1) != sha(tifffile.imread(c["gray_crop_path"])):
                raise RuntimeError(f"f=1 resample differs from dev crop {c['crop_id']}")
            for f in FACTORS:
                cid = f"{c['crop_id']}_{tag(f)}"
                d = ROOT / tag(f) / cid
                d.mkdir(parents=True, exist_ok=True)
                g, box = resample(gray_v, centre, f, 1, gray_v.dtype)
                lab, _ = resample(gt_v, centre, f, 0, gt_v.dtype)
                tifffile.imwrite(d / "gray.tif", g)
                tifffile.imwrite(d / "labels.tif", lab)
                mgv = int(round(1000 * f ** 3))
                prot = reference_protocol(lab, mgv, True)
                n_valid = int(len(prot["valid_ids"]))
                crops_out.append({"crop_id": cid, "dev_crop_id": c["crop_id"], "scan": scan, "specimen": "EFR",
                                  "factor": f, "centre_zyx": centre, "source_box_zyx": box, "size": SIZE,
                                  "gray_crop_path": str(d / "gray.tif"), "gt_path": str(d / "labels.tif"),
                                  "min_gt_voxels": mgv})
                rows.append({"crop_id": cid, "factor": f, "dev_crop_id": c["crop_id"], "gray_sha256": sha(g),
                             "labels_sha256": sha(lab), "shape": "x".join(map(str, g.shape)), "gray_dtype": str(g.dtype),
                             "labels_dtype": str(lab.dtype), "n_labels": int(len(np.unique(lab)) - (lab == 0).any()),
                             "min_gt_voxels": mgv, "valid_gt": n_valid})
                print(f"[synth] {cid}: labels {rows[-1]['n_labels']}, valid {n_valid} (min {mgv} vox)", flush=True)
    meta = {"schema_version": 1, "created": "2026-11-07", "description": __doc__, "factors": list(FACTORS),
            "source": str(DEV), "crop_convention": "zyx, vol[z:z+200, y:y+200, x:x+200]", "crops": crops_out}
    OUT_JSON.write_text(json.dumps(meta, indent=1), encoding="utf-8")
    with OUT_CSV.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"[done] {len(crops_out)} synthetic crops -> {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
