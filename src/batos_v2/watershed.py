"""Stage 4: fixed-mask marker-controlled distance watershed + prediction size filter variants.

Recovered and ported exactly:
- baseline_existing markers = batos_markers.py with the run_batos_pipeline.sh flags (used ONLY to
  enlarge the frozen mask, as in E11/E23c: --mask-marker ToS MIN + --mask-marker baseline)
- fixed mask + watershed = e11_run_fixed_mask_watershed.py: Otsu (no offset) on the p1-p99
  normalized openrec image, 6-conn closing x1, fill holes, components >= 500, OR mask-marker voxels;
  relief = -(EDT normalized p1-p99 inside the mask); markers clipped to the mask and relabeled
  (6-conn); skimage watershed, connectivity 1, no watershed line.
Reconstructed: BG10 size filter (labels < N voxels removed; surviving labels keep their ids so each
basin stays linked to its seed for seed-ignore evaluation).

Outputs: frozen_mask_uint8.tif, negative_edt_float32.tif, baseline_markers.tif,
normalized_markers/<label>.tif, segmentations/<label>_distance_c1[_minN].tif, watershed_results.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.segmentation import watershed

from .common import CONNECTIVITY_6, array_sha256, load_config, normalize_p1_p99, read3d, robust_u8, write3d, write_json
from .gate import clean_bright_otsu_mask


# ============================================================ baseline markers (batos_markers.py)
def _offsets_sphere(r):
    r = int(r)
    return [(dz, dy, dx) for dz in range(-r, r + 1) for dy in range(-r, r + 1) for dx in range(-r, r + 1)
            if dz * dz + dy * dy + dx * dx <= r * r]


def _offsets_ring2d(inner, outer):
    i2, o2, r = int(inner) ** 2, int(outer) ** 2, int(outer)
    return [(dy, dx) for dy in range(-r, r + 1) for dx in range(-r, r + 1) if i2 < dx * dx + dy * dy <= o2]


def _sample_sphere(vol, z, y, x, offsets):
    d, h, w = vol.shape
    vals = [vol[z + dz, y + dy, x + dx] for dz, dy, dx in offsets
            if 0 <= z + dz < d and 0 <= y + dy < h and 0 <= x + dx < w]
    return np.asarray(vals)


def _sample_ring2d(vol, z, y, x, offsets):
    h, w = vol.shape[1:]
    plane = vol[z]
    vals = [plane[y + dy, x + dx] for dy, dx in offsets if 0 <= y + dy < h and 0 <= x + dx < w]
    return np.asarray(vals)


def baseline_markers(flat_u8: np.ndarray, internal: np.ndarray, cfg: dict[str, Any]) -> tuple[np.ndarray, dict]:
    structure = ndi.generate_binary_structure(3, cfg["connectivity"])
    lab, n = ndi.label(internal, structure=structure)
    slices = ndi.find_objects(lab)
    center_sphere = _offsets_sphere(cfg["center_radius"])
    ring2d = _offsets_ring2d(cfg["ring_inner"], cfg["ring_outer"])
    D, H, W = internal.shape
    cands = []
    for cid, slc in enumerate(slices, start=1):
        if slc is None:
            continue
        comp = lab[slc] == cid
        area = int(comp.sum())
        if area < cfg["area_min"] or area > cfg["area_max"]:
            continue
        dz, dy, dx = (s.stop - s.start for s in slc)
        max_dim = max(dx, dy, dz)
        min_dim = max(1, min(dx, dy, dz))
        elong = max_dim / min_dim
        fill = area / max(1, dx * dy * dz)
        if max_dim > cfg["max_dim"] or elong > cfg["max_elongation"] or fill < cfg["min_fill"]:
            continue
        coords = np.argwhere(comp)
        cz, cy, cx = np.rint(coords.mean(axis=0) + np.asarray([s.start for s in slc])).astype(int)
        cz, cy, cx = int(np.clip(cz, 0, D - 1)), int(np.clip(cy, 0, H - 1)), int(np.clip(cx, 0, W - 1))
        center_vals = _sample_sphere(flat_u8, cz, cy, cx, center_sphere)
        ring_vals = _sample_ring2d(flat_u8, cz, cy, cx, ring2d)
        center_mean = float(center_vals.mean()) if center_vals.size else 255.0
        ring_mean = float(ring_vals.mean()) if ring_vals.size else 0.0
        ring_contrast = ring_mean - center_mean
        if center_mean > cfg["center_dark_max"] or ring_contrast < cfg["min_ring_contrast"]:
            continue
        score = ring_contrast + max(0.0, 205.0 - center_mean) * 0.02
        cands.append({"slice": slc, "mask": comp, "area": area, "elongation": elong,
                      "z": cz, "y": cy, "x": cx, "score": score})
    # nms() of batos_markers.py
    cands = sorted(cands, key=lambda c: (-c["score"], c["area"], c["elongation"]))
    accepted = []
    md2 = float(cfg["nms_min_dist"] ** 2)
    for c in cands:
        p = np.asarray([c["z"], c["y"], c["x"]], dtype=np.float64)
        if all(float(np.sum((p - np.asarray([a["z"], a["y"], a["x"]], dtype=np.float64)) ** 2)) >= md2
               for a in accepted):
            accepted.append(c)
            if cfg["max_markers"] > 0 and len(accepted) >= cfg["max_markers"]:
                break
    out = np.zeros(internal.shape, dtype=np.uint8)
    for c in accepted:
        out[c["slice"]][c["mask"]] = 255
    return out, {"internal_components": int(n), "candidates_before_nms": len(cands), "kept": len(accepted)}


# ============================================================ fixed mask + watershed (e11)
def normalize01_in_mask(array: np.ndarray, mask: np.ndarray) -> np.ndarray:
    values = array[mask & np.isfinite(array)]
    if values.size == 0:
        return np.zeros(array.shape, dtype=np.float32)
    low, high = np.percentile(values, [1.0, 99.0])
    if high <= low:
        low, high = float(values.min()), float(values.max())
    if high <= low:
        return np.zeros(array.shape, dtype=np.float32)
    result = np.clip((array - low) / (high - low), 0.0, 1.0)
    result[~np.isfinite(result)] = 0.0
    return result.astype(np.float32, copy=False)


def frozen_mask(flat_image: np.ndarray, mask_markers: list[np.ndarray], cfg: dict[str, Any]):
    mask, threshold = clean_bright_otsu_mask(normalize_p1_p99(flat_image), cfg["mask_min_component"])
    for extra in mask_markers:
        mask |= extra > 0
    if not mask.any():
        raise RuntimeError("Frozen mask is empty.")
    return mask, threshold


def normalize_markers_in_mask(raw: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, int]:
    positive = raw > 0
    outside = int((positive & ~mask).sum())
    positive &= mask
    labels, count = ndi.label(positive, structure=CONNECTIVITY_6)
    if count == 0:
        raise ValueError("Marker volume empty after applying the mask.")
    return labels.astype(np.int32, copy=False), outside


def size_filter(labels: np.ndarray, min_voxels: int) -> np.ndarray:
    counts = np.bincount(labels.ravel())
    remove = counts < int(min_voxels)
    remove[0] = False
    out = labels.copy()
    out[remove[labels]] = 0
    return out


def run(flat_image: np.ndarray, selector_u8: np.ndarray, internal: np.ndarray, marker_paths: dict[str, Path],
        config: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    out_dir = Path(out_dir)
    seg_dir, norm_dir = out_dir / "segmentations", out_dir / "normalized_markers"
    seg_dir.mkdir(parents=True, exist_ok=True)
    norm_dir.mkdir(parents=True, exist_ok=True)
    cfg = config["watershed"]

    baseline, baseline_info = baseline_markers(robust_u8(selector_u8), internal > 0, config["baseline_markers"])
    write3d(out_dir / "baseline_markers.tif", baseline)
    union_sources = {"tos_min": read3d(marker_paths["tos_e11_min"]), "baseline_existing": baseline}
    mask, threshold = frozen_mask(flat_image, [union_sources[k] for k in cfg["union_markers"]], cfg)
    distance = ndi.distance_transform_edt(mask)
    surface = -normalize01_in_mask(distance, mask)
    write3d(out_dir / "frozen_mask_uint8.tif", mask.astype(np.uint8))
    write3d(out_dir / "negative_edt_float32.tif", surface.astype(np.float32))

    rows = []
    for label, path in marker_paths.items():
        markers, outside = normalize_markers_in_mask(read3d(path, f"markers {label}"), mask)
        norm_path = write3d(norm_dir / f"{label}.tif", markers)
        prediction = watershed(surface, markers=markers, connectivity=cfg["connectivity"], mask=mask,
                               watershed_line=False).astype(np.int32, copy=False)
        for mode in config["size_filter"]["modes"]:
            mode = "off" if mode is False or str(mode).lower() == "off" else int(mode)
            if mode == "off":
                seg, suffix = prediction, ""
            else:
                seg, suffix = size_filter(prediction, int(mode)), f"_min{int(mode)}"
            seg_path = write3d(seg_dir / f"{label}_distance_c1{suffix}.tif", seg)
            rows.append({"configuration": f"{label}_distance_c1{suffix}", "marker_source": label,
                         "size_filter": str(mode), "marker_count": int(markers.max()),
                         "marker_voxels_outside_frozen_mask": outside,
                         "prediction_count": int(len(np.unique(seg)) - (1 if (seg == 0).any() else 0)),
                         "normalized_marker": str(norm_path), "segmentation": str(seg_path),
                         "segmentation_array_sha256": array_sha256(seg)})
    frame = pd.DataFrame(rows)
    frame.to_csv(out_dir / "watershed_results.csv", index=False)
    info = {"otsu_threshold": threshold, "mask_fraction": float(mask.mean()),
            "mask_array_sha256": array_sha256(mask.astype(np.uint8)), "baseline_markers": baseline_info,
            "results_csv": str(out_dir / "watershed_results.csv"), "rows": rows}
    write_json(out_dir / "watershed_manifest.json", info)
    return info


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="batos-paper-v2 stage 4: fixed-mask distance watershed")
    parser.add_argument("--flat-image", type=Path, required=True, help="*_flatten_openrec_r*_8u.tif")
    parser.add_argument("--selector", type=Path, required=True, help="*_flatten_lambda1p0_8u.tif")
    parser.add_argument("--internal-dark", type=Path, required=True)
    parser.add_argument("--marker", action="append", required=True, metavar="LABEL=TIFF",
                        help="must include tos_e11_min=...")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("configs/batos_paper_v2.yaml"))
    args = parser.parse_args(argv)
    markers = dict(item.split("=", 1) for item in args.marker)
    if "tos_e11_min" not in markers:
        raise ValueError("--marker tos_e11_min=... is required (it defines the frozen mask)")
    info = run(read3d(args.flat_image), read3d(args.selector), read3d(args.internal_dark),
               {k: Path(v) for k, v in markers.items()}, load_config(args.config), args.out_dir)
    print(pd.DataFrame(info["rows"])[["configuration", "marker_count", "prediction_count"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
