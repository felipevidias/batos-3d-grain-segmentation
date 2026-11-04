"""Stage 3: support-collapse gate and selective marker routing (BA-TOS_sel).

Recovered and ported exactly:
- gate statistics = e23b/e23c compute_gate: merge markers.csv (node) with
  candidate_diagnostics_no_gt.csv one-to-one, take region_voxels of the selected marker nodes,
  g = median <= 3 and q75 <= 4 (pandas linear quantile)
- fallback markers = e11_build_control_markers.py --mode otsu_edt --nms-radius 6
  --match-count-from <ToS markers> on the openrec ("FLAT") image
- routing = e23c: gate true -> EDT-Otsu, else ToS MIN

Outputs: markers_edt_otsu.tif / .csv, batos_selective_markers_uint16.tif, route.json
"""

from __future__ import annotations

import argparse
import itertools
import math
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu

from .common import CONNECTIVITY_6, array_sha256, load_config, normalize_p1_p99, read3d, write3d, write_json


# ============================================================ gate (e23b/e23c)
def compute_gate(diagnostics: pd.DataFrame, markers: pd.DataFrame, cfg: dict[str, Any]) -> dict[str, Any]:
    missing = {"node", "polarity", "score", "has_reconstructed_region", "region_voxels"} - set(diagnostics.columns)
    if missing:
        raise ValueError(f"diagnostics missing columns: {sorted(missing)}")
    diagnostics = diagnostics.copy()
    diagnostics["node"] = pd.to_numeric(diagnostics["node"], errors="raise").astype(np.int64)
    marker_nodes = markers[["node"]].copy()
    marker_nodes["node"] = pd.to_numeric(marker_nodes["node"], errors="raise").astype(np.int64)
    selected = marker_nodes.merge(diagnostics, on="node", how="left", validate="one_to_one")
    if selected["region_voxels"].isna().any():
        raise RuntimeError("Selected marker nodes missing from diagnostics.")
    values = pd.to_numeric(selected["region_voxels"], errors="coerce").dropna()
    if values.empty:
        raise RuntimeError("No valid region_voxels for the gate.")
    median_value = float(values.median())
    q25_value = float(values.quantile(0.25))
    q75_value = float(values.quantile(0.75))
    return {
        "gate_name": "region_median_le3_and_q75_le4",
        "gate_true": bool(median_value <= cfg["median_max"] and q75_value <= cfg["q75_max"]),
        "region_voxels_q25": q25_value,
        "region_voxels_median": median_value,
        "region_voxels_q75": q75_value,
        "median_limit": float(cfg["median_max"]),
        "q75_limit": float(cfg["q75_max"]),
        "n_markers": int(len(values)),
        "gt_used": False,
    }


# ============================================================ EDT-Otsu control markers (e11_build_control_markers)
def clean_bright_otsu_mask(image_norm: np.ndarray, min_component: int = 500) -> tuple[np.ndarray, float]:
    threshold = float(threshold_otsu(image_norm))
    mask = image_norm > threshold
    mask = ndi.binary_closing(mask, structure=CONNECTIVITY_6, iterations=1)
    mask = ndi.binary_fill_holes(mask)
    labels, count = ndi.label(mask, structure=CONNECTIVITY_6)
    if count == 0:
        raise RuntimeError("Otsu bright mask is empty.")
    sizes = np.bincount(labels.ravel())
    keep = sizes >= int(min_component)
    keep[0] = False
    return keep[labels], threshold


def _greedy_nms_points(frame: pd.DataFrame, radius: float) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    points = frame[["z", "y", "x"]].to_numpy(dtype=np.int64)
    physical = points.astype(float)
    priority = frame["priority"].to_numpy(dtype=float)
    source_index = frame["source_index"].to_numpy(dtype=np.int64)
    order = np.lexsort((points[:, 2], points[:, 1], points[:, 0], source_index, -priority))
    radius_squared = float(radius * radius)
    grid: dict[tuple[int, int, int], list[int]] = {}
    kept: list[int] = []
    offsets = list(itertools.product((-1, 0, 1), repeat=3))
    for index in order:
        point = physical[index]
        cell = tuple(int(math.floor(v / radius)) for v in point)
        reject = False
        for off in offsets:
            for previous in grid.get((cell[0] + off[0], cell[1] + off[1], cell[2] + off[2]), []):
                delta = point - physical[previous]
                if float(np.dot(delta, delta)) < radius_squared:
                    reject = True
                    break
            if reject:
                break
        if reject:
            continue
        kept.append(int(index))
        grid.setdefault(cell, []).append(int(index))
    return frame.iloc[kept].copy().reset_index(drop=True)


def marker_component_count(volume: np.ndarray) -> int:
    _, count = ndi.label(volume > 0, structure=CONNECTIVITY_6)
    if count <= 0:
        raise ValueError("Reference marker volume is empty.")
    return int(count)


def edt_otsu_markers(image: np.ndarray, nms_radius: float, target_count: int | None):
    mask, threshold = clean_bright_otsu_mask(normalize_p1_p99(image))
    edt = ndi.distance_transform_edt(mask)
    coords = peak_local_max(edt, min_distance=1, threshold_abs=0.0, exclude_border=False,
                            num_peaks=np.inf, labels=mask.astype(np.uint8), p_norm=2)
    if coords.size == 0:
        raise RuntimeError("No EDT local maximum found.")
    values = edt[tuple(coords[:, a] for a in range(3))]
    candidates = pd.DataFrame({"source_index": np.arange(len(coords), dtype=np.int64),
                               "z": coords[:, 0], "y": coords[:, 1], "x": coords[:, 2],
                               "priority": values, "control_mode": "otsu_edt"})
    before = len(candidates)
    selected = _greedy_nms_points(candidates, nms_radius)
    after = len(selected)
    selected = selected.sort_values(["priority", "z", "y", "x", "source_index"],
                                    ascending=[False, True, True, True, True]).reset_index(drop=True)
    if target_count is not None:
        if len(selected) < target_count:
            raise RuntimeError(f"otsu_edt: only {len(selected)} markers after NMS; target={target_count}")
        selected = selected.iloc[:target_count].copy().reset_index(drop=True)
    selected["marker_id"] = np.arange(1, len(selected) + 1, dtype=np.int64)
    volume = np.zeros(image.shape, dtype=np.uint16 if len(selected) <= 65535 else np.uint32)
    for label, point in enumerate(selected[["z", "y", "x"]].to_numpy(np.int64), 1):
        volume[tuple(point)] = label
    info = {"otsu_threshold": threshold, "otsu_mask_fraction": float(mask.mean()),
            "candidate_count_before_nms": before, "marker_count_after_nms": after,
            "target_count": target_count, "final_marker_count": int(len(selected))}
    return selected, volume, info


# ============================================================ stage runner
def run(diagnostics_csv: Path, markers_csv: Path, tos_markers_tif: Path, flat_image: np.ndarray,
        config: dict[str, Any], out_dir: Path, flat_image_path: Path | None = None) -> dict[str, Any]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gate = compute_gate(pd.read_csv(diagnostics_csv), pd.read_csv(markers_csv), config["gate"])
    tos_volume = read3d(tos_markers_tif, "ToS markers")
    target = marker_component_count(tos_volume)
    edt_frame, edt_volume, edt_info = edt_otsu_markers(flat_image, config["gate"]["edt_nms_radius"], target)
    edt_frame.to_csv(out_dir / "markers_edt_otsu_port.csv", index=False)
    edt_tif = write3d(out_dir / "markers_edt_otsu_port.tif", edt_volume)
    if config["gate"].get("edt_implementation", "port") == "e11_verbatim":
        if flat_image_path is None:
            raise ValueError("gate.edt_implementation=e11_verbatim needs flat_image_path")
        from .recovered import e11_edt_otsu
        verbatim = e11_edt_otsu(config, Path(flat_image_path), Path(tos_markers_tif), out_dir / "markers_edt_otsu_e11")
        verbatim_volume = read3d(verbatim["markers_tif"])
        edt_info = dict(edt_info, **verbatim, port_identical_to_verbatim=bool(np.array_equal(verbatim_volume, edt_volume)))
        edt_tif, edt_volume = Path(verbatim["markers_tif"]), verbatim_volume
    source, source_path = ("edt_otsu", edt_tif) if gate["gate_true"] else ("tos_e11_min", Path(tos_markers_tif))
    selective = out_dir / "batos_selective_markers_uint16.tif"
    shutil.copy2(source_path, selective)
    route = {"gate": gate, "rule": {"gate_false": "tos_e11_min", "gate_true": "edt_otsu"},
             "selected_source": source, "selected_source_path": str(source_path),
             "selective_markers": str(selective), "edt_otsu": edt_info,
             "edt_otsu_markers": str(edt_tif), "tos_markers": str(tos_markers_tif),
             "marker_counts": {"tos_e11_min": target, "edt_otsu": int(edt_volume.max())},
             "array_sha256": {"tos_e11_min": array_sha256(tos_volume), "edt_otsu": array_sha256(edt_volume)}}
    write_json(out_dir / "route.json", route)
    return route


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="batos-paper-v2 stage 3: support-collapse gate and routing")
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--markers-csv", type=Path, required=True)
    parser.add_argument("--markers-tif", type=Path, required=True)
    parser.add_argument("--flat-image", type=Path, required=True, help="*_flatten_openrec_r*_8u.tif")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("configs/batos_paper_v2.yaml"))
    args = parser.parse_args(argv)
    route = run(args.diagnostics, args.markers_csv, args.markers_tif, read3d(args.flat_image),
                load_config(args.config), args.out_dir)
    g = route["gate"]
    print(f"[gate] median={g['region_voxels_median']} q75={g['region_voxels_q75']} "
          f"gate={g['gate_true']} -> {route['selected_source']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
