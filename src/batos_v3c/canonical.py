"""BA-TOS v3c: canonical-scale normalization around the frozen v2-fixflat pipeline (Phase 5i, docs/V3C_CANONICAL.md).

r_granulo (batos_v3.scale, native crop) -> f = R_REF / r_granulo -> cubic resample of the grey crop to canonical scale
-> batos_v3.pipeline.run_crop(preset="fixflat") with e13 replaced by (r = R_REF, openrec = max(8, round(1.5 R_REF)))
-> size filter 500 at canonical scale -> nearest-neighbour map back of labels, coordinate map back of markers.
GT is never read. r_granulo is the only scale input.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import time
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

from batos_v2.run_crop import git_state
from batos_v3.constants import R_REF
from batos_v3.pipeline import SEGMENTATIONS, StageError, run_crop
from batos_v3.scale import estimate

SIZE_FILTER = 500  # canonical voxels (frozen v2 value)
ROUTES = list(SEGMENTATIONS)


def canonical_shape(shape, f: float) -> tuple[int, ...]:
    return tuple(max(1, int(round(n * f))) for n in shape)


def resample_in(gray: np.ndarray, cshape) -> np.ndarray:
    zoom = [c / n for c, n in zip(cshape, gray.shape)]
    out = ndi.zoom(gray.astype(np.float32), zoom, output=np.float32, order=3, mode="grid-mirror", grid_mode=True)
    if out.shape != tuple(cshape):
        raise RuntimeError(f"zoom shape {out.shape} != {cshape}")
    if np.issubdtype(gray.dtype, np.integer):
        info = np.iinfo(gray.dtype)
        out = np.clip(np.rint(out), info.min, info.max)
    return out.astype(gray.dtype)


def labels_back(lab: np.ndarray, nshape) -> np.ndarray:
    """Nearest neighbour: native voxel i <- canonical voxel floor((i + 0.5) * Nc / N)."""
    idx = [np.minimum(((np.arange(n) + 0.5) * (c / n)).astype(np.int64), c - 1) for n, c in zip(nshape, lab.shape)]
    return lab[np.ix_(*idx)]


def markers_back(mk: np.ndarray, nshape) -> np.ndarray:
    """Each canonical marker voxel j -> native round((j + 0.5) * N / Nc - 0.5), label kept."""
    out = np.zeros(nshape, mk.dtype)
    pos = np.nonzero(mk)
    if not len(pos[0]):
        return out
    p = [np.clip(np.rint((j + 0.5) * (n / c) - 0.5).astype(np.int64), 0, n - 1)
         for j, n, c in zip(pos, nshape, mk.shape)]
    out[tuple(p)] = mk[pos]
    return out


def size_filter(pred: np.ndarray, min_size: int) -> np.ndarray:
    ids, counts = np.unique(pred, return_counts=True)
    keep = ids[(ids > 0) & (counts >= min_size)]
    return np.where(np.isin(pred, keep), pred, 0).astype(pred.dtype)


def run(crop_id: str, gray_path: Path, out: Path, r_probe: float | None = None) -> dict:
    """Writes out/{scale.json, canon/..., native/<route>_{labels,labels_own,markers}.tif, v3c_run.json, DONE|FAILED.json}.

    r_probe: r_granulo already estimated on the 200^3 probe that defined this window (canonical field-of-view protocol);
    when None, r_granulo is estimated on the crop itself (a 200^3 crop is its own probe)."""
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    gray = tifffile.imread(gray_path)
    if r_probe is None:
        sc, _ = estimate(gray)
        sc = {k: v for k, v in sc.items() if not isinstance(v, dict)} | {"source": "estimated on this crop"}
    else:
        sc = {"r_granulo": float(r_probe), "source": "200^3 probe at the same centre (window defined from it)"}
    r = float(sc["r_granulo"])
    f = R_REF / r
    (out / "scale.json").write_text(json.dumps(sc, indent=1), encoding="utf-8")
    cshape = canonical_shape(gray.shape, f)
    rec = {"crop_id": crop_id, "method": "v3c", "gray": str(gray_path), "git": git_state(),
           "timestamp": dt.datetime.now().astimezone().isoformat(timespec="seconds"), "r_granulo": r, "f": f,
           "native_shape": list(gray.shape), "canonical_shape": list(cshape), "size_filter_canonical": SIZE_FILTER,
           "scale_source": sc["source"], "gt_used": False}
    cdir = out / "canon"
    cdir.mkdir(exist_ok=True)
    if cshape == gray.shape:
        cgray = Path(gray_path)  # f = 1 (or rounds to the same grid): no interpolation
        rec["resampled"] = False
    else:
        cgray = cdir / "gray_canonical.tif"
        tifffile.imwrite(cgray, resample_in(gray, cshape), compression="zlib")
        rec["resampled"] = True
    del gray
    params = cdir / "canonical_params.json"
    params.write_text(json.dumps({"estimated_grain_radius": R_REF, "openrec_radius": int(max(8, round(1.5 * R_REF))),
                                  "source": "v3c canonical scale (replaces e13; r_granulo is the only scale input)"}),
                      encoding="utf-8")
    try:
        crec = run_crop(crop_id, cgray, cdir / "run", "fixflat", params)
    except StageError as e:
        rec |= {"error": str(e), "seconds_total": round(time.time() - t0, 1)}
        (out / "FAILED.json").write_text(json.dumps(rec, indent=1, default=str), encoding="utf-8")
        return rec
    rec |= {"canonical_run": {k: crec.get(k) for k in ("openrec_radius", "bg_ref", "gate", "tos_markers",
                                                       "seconds_total", "peak_memory_mb")}}
    ndir = out / "native"
    ndir.mkdir(exist_ok=True)
    for route, (ws, label) in SEGMENTATIONS.items():
        seg = tifffile.imread(cdir / "run" / ws / "segmentations" / f"{label}_distance_c1.tif")
        mk = tifffile.imread(cdir / "run" / ws / "normalized_markers" / f"{label}.tif")
        tifffile.imwrite(ndir / f"{route}_labels.tif", labels_back(seg, rec["native_shape"]), compression="zlib")
        tifffile.imwrite(ndir / f"{route}_labels_own.tif", labels_back(size_filter(seg, SIZE_FILTER),
                                                                       rec["native_shape"]), compression="zlib")
        tifffile.imwrite(ndir / f"{route}_markers.tif", markers_back(mk, rec["native_shape"]), compression="zlib")
    rec["seconds_total"] = round(time.time() - t0, 1)
    rec["peak_memory_mb"] = crec.get("peak_memory_mb")
    (out / "v3c_run.json").write_text(json.dumps(rec, indent=1, default=str), encoding="utf-8")
    (out / "DONE").write_text("ok", encoding="utf-8")
    return rec


def outputs(out: Path, shape) -> dict:
    """route -> (labels none, labels own, markers) in native space; zeros for a failed run."""
    if (out / "FAILED.json").exists():
        z = np.zeros(shape, np.uint32)
        return {r: (z, z, z) for r in ROUTES}
    n = out / "native"
    return {r: tuple(tifffile.imread(n / f"{r}_{k}.tif") for k in ("labels", "labels_own", "markers")) for r in ROUTES}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--crop-id", required=True)
    ap.add_argument("--gray", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    rec = run(a.crop_id, a.gray, a.out)
    print(json.dumps({k: rec.get(k) for k in ("crop_id", "r_granulo", "f", "canonical_shape", "error")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
