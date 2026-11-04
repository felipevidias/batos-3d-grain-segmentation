"""Stage 2: boundary rejection, ToS candidate selection, spatial diagnostics, internal markers.

Recovered and ported exactly:
- internal_dark mask = batos_background.py (choose_dark_threshold + border propagation): the
  operational boundary rule (dark components connected to any of the 6 faces are external).
- candidate scoring = batos_tos_seed.py (origin/batos-ls-auto-tile-pilot), Higra 3D ToS.
- marker consolidation + NMS = e11_build_markers_ablation.py.
Reconstructed (lost scripts batos_tos_seed_true_dual.py, e10_spatialize_candidates_no_gt.py):
- BG8 true-dual split (top K/2 MIN-like, top K/2 MAX-like; MAX scoring mirrors MIN scoring)
- BG9 exclusive supports (each leaf voxel -> nearest selected ancestor), BG12 representative point

Outputs in --out-dir:
  <prefix>_internal_dark_candidates.tif           uint8 0/255
  <prefix>_tosconop_autodual_candidates.csv       selected nodes (both polarities)
  candidate_diagnostics_no_gt.csv                 columns read by E11/E23 scripts
  markers.csv, markers_uint16.tif                 ToS MIN point markers (labels 1..N)
  tos_select.json                                 statistics
"""

from __future__ import annotations

import argparse
import itertools
import math
import time
from pathlib import Path
from typing import Any

import higra as hg
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from .common import CONNECTIVITY_6, load_config, read3d, robust_u8, write3d, write_json


# ============================================================ boundary rejection (batos_background.py)
def otsu_u8(vol: np.ndarray) -> int:
    hist = np.bincount(vol.ravel(), minlength=256).astype(np.float64)
    total = hist.sum()
    if total <= 0:
        return 0
    p = hist / total
    omega = np.cumsum(p)
    mu = np.cumsum(p * np.arange(256))
    mt = mu[-1]
    denom = omega * (1.0 - omega)
    denom[denom <= 1e-12] = np.nan
    sigma = (mt * omega - mu) ** 2 / denom
    return int(np.nanargmax(sigma))


def internal_dark_mask(flat_u8: np.ndarray, cfg: dict[str, Any]) -> tuple[np.ndarray, dict[str, Any]]:
    nz = flat_u8[flat_u8 > 0]
    if nz.size == 0:
        nz = flat_u8.ravel()
    pval = float(np.percentile(nz, cfg["dark_percentile"]))
    otsu = otsu_u8(flat_u8)
    thr = int(np.clip(min(otsu + cfg["otsu_offset"], pval, cfg["dark_threshold_max"]), 1, 255))
    dark = flat_u8 <= thr
    border = np.zeros_like(dark)
    for index in [(0,), (-1,)]:
        border[index] |= dark[index]
    border[:, 0, :] |= dark[:, 0, :]
    border[:, -1, :] |= dark[:, -1, :]
    border[:, :, 0] |= dark[:, :, 0]
    border[:, :, -1] |= dark[:, :, -1]
    structure = ndi.generate_binary_structure(3, cfg["connectivity"])
    external = ndi.binary_propagation(border, structure=structure, mask=dark)
    internal = dark & ~external
    return internal, {"dark_threshold": thr, "otsu_u8": otsu, "dark_percentile_value": pval,
                      "dark_voxels": int(dark.sum()), "external_dark_voxels": int(external.sum()),
                      "internal_dark_voxels": int(internal.sum())}


# ============================================================ ToS candidates (batos_tos_seed.py + BG8)
def _clamp01(x):
    return np.clip(x, 0.0, 1.0)


def _nearest_support(node, parents, altitudes, area, root, polarity_fn, cfg):
    """Port of nearest_max_support (and its mirror for MAX candidates)."""
    cur = int(parents[node])
    steps = 1
    while cur != root and cur >= 0 and steps <= cfg["support_max_steps"]:
        if polarity_fn(cur) and cfg["support_area_min"] <= int(area[cur]) <= cfg["support_area_max"]:
            return cur, steps
        nxt = int(parents[cur])
        if nxt == cur:
            break
        cur = nxt
        steps += 1
    return -1, steps


def score_candidates(image_u8: np.ndarray, cfg: dict[str, Any], log=print):
    t0 = time.time()
    tree, altitudes = hg.component_tree_tree_of_shapes_image3d(image_u8)
    altitudes = np.asarray(altitudes).astype(np.int64)
    parents = tree.parents().astype(np.int64)
    root = int(tree.root())
    area = hg.attribute_area(tree).astype(np.int64)
    height = hg.attribute_height(tree, altitudes).astype(np.float64)
    n_nodes = len(parents)
    log(f"[tos] tree nodes={n_nodes} leaves={tree.num_leaves()} built in {time.time() - t0:.1f}s")

    not_root = np.arange(n_nodes) != root
    parent_alt = altitudes[parents]
    is_min = not_root & (altitudes < parent_alt)
    is_max = not_root & (altitudes > parent_alt)
    min_like = np.flatnonzero(is_min)
    max_like = np.flatnonzero(is_max)
    max_min_height = float(height[min_like].max()) if min_like.size else 1.0
    max_max_height = float(height[max_like].max()) if max_like.size else 1.0
    avg_min_area = float(area[min_like].mean()) if min_like.size else 1.0
    avg_max_area = float(area[max_like].mean()) if max_like.size else 1.0
    nz = image_u8[image_u8 > 0]
    dark_alt_max = float(np.percentile(nz, cfg["dark_percentile_max"])) if nz.size else 255.0
    bright_alt_min = float(np.percentile(nz, 100.0 - cfg["dark_percentile_max"])) if nz.size else 0.0

    log_min = math.log(float(max(1, cfg["area_min"])) + 1.0)
    log_max = math.log(float(max(cfg["area_min"] + 1, cfg["area_max"])) + 1.0)

    def score_polarity(polarity: str):
        if polarity == "min":
            nodes, avg_area, max_h, max_h_other = min_like, avg_min_area, max_min_height, max_max_height
            delta = parent_alt - altitudes
            extreme_ok = altitudes <= dark_alt_max
            extreme_score = _clamp01((dark_alt_max - altitudes) / max(1.0, dark_alt_max))
            support_fn = lambda n: bool(is_max[n])
        else:
            nodes, avg_area, max_h, max_h_other = max_like, avg_max_area, max_max_height, max_min_height
            delta = altitudes - parent_alt
            extreme_ok = altitudes >= bright_alt_min
            extreme_score = _clamp01((altitudes - bright_alt_min) / max(1.0, 255.0 - bright_alt_min))
            support_fn = lambda n: bool(is_min[n])
        a = area[nodes]
        keep = ((a >= cfg["area_min"]) & (a <= cfg["area_max"]) & (a <= cfg["area_factor"] * avg_area)
                & extreme_ok[nodes] & (delta[nodes] >= cfg["delta_parent_min"]))
        nodes = nodes[keep]
        h_score = _clamp01(height[nodes] / max(1e-9, cfg["height_factor"] * max_h))
        d_score = _clamp01(delta[nodes] / max(1.0, cfg["delta_norm"]))
        area_pos = (np.log(area[nodes].astype(np.float64) + 1.0) - log_min) / max(1e-9, log_max - log_min)
        a_score = _clamp01(1.0 - np.abs(area_pos - cfg["area_preferred_position"]) / max(1e-9, cfg["area_preferred_width"]))
        x_score = extreme_score[nodes]
        partial = (cfg["w_height"] * h_score + cfg["w_delta"] * d_score + cfg["w_area"] * a_score
                   + cfg["w_dark"] * x_score)
        # support_score <= 1: nodes that cannot reach min_score even with full support are skipped (same result)
        reachable = partial + cfg["w_support"] >= cfg["min_score"]
        rows = []
        for i in np.flatnonzero(reachable):
            node = int(nodes[i])
            support, steps = _nearest_support(node, parents, altitudes, area, root, support_fn, cfg)
            s_score = 0.0
            s_area = 0
            s_height = 0.0
            s_contrast = 0.0
            if support >= 0:
                s_area = int(area[support])
                s_height = float(height[support])
                s_contrast = abs(float(altitudes[support]) - float(altitudes[node]))
                sh = min(1.0, max(0.0, s_height / max(1e-9, 0.10 * max_h_other)))
                ratio = float(s_area) / max(1.0, float(area[node]))
                rs = min(1.0, max(0.0, 1.0 - abs(math.log(ratio + 1e-9) - math.log(cfg["support_area_ratio_preferred"]))
                                  / max(1e-9, math.log(cfg["support_area_ratio_width"]))))
                cs = min(1.0, max(0.0, s_contrast / max(1.0, cfg["support_contrast_norm"])))
                ds = min(1.0, max(0.0, 1.0 - (steps - 1) / max(1.0, cfg["support_max_steps"])))
                s_score = 0.35 * sh + 0.30 * rs + 0.25 * cs + 0.10 * ds
            score = float(partial[i] + cfg["w_support"] * s_score)
            if score < cfg["min_score"]:
                continue
            rows.append({
                "node": node, "polarity": polarity, "score": score, "area": int(area[node]),
                "height": float(height[node]), "alt": int(altitudes[node]),
                "parent_alt": int(parent_alt[node]), "delta_parent": int(delta[node]),
                "support": int(support), "support_area": s_area, "support_height": s_height,
                "support_contrast": s_contrast, "support_steps": int(steps),
                "h_score": float(h_score[i]), "delta_score": float(d_score[i]),
                "area_score": float(a_score[i]), "dark_score": float(x_score[i]),
                "support_score": s_score,
            })
        # batos_tos_seed.py ordering: sort by (score, height, -area) descending, stable on node order
        rows.sort(key=lambda r: (r["score"], r["height"], -r["area"]), reverse=True)
        return rows

    t1 = time.time()
    min_rows = score_polarity("min")
    max_rows = score_polarity("max")
    k_total = int(cfg["k_total"])
    if cfg["k_split"] == "half":
        selected = min_rows[: k_total // 2] + max_rows[: k_total - k_total // 2]
    elif cfg["k_split"] == "min_only":
        selected = min_rows[:k_total]
    else:
        raise ValueError(f"Unknown tos.k_split {cfg['k_split']!r}")
    info = {
        "num_tree_nodes": n_nodes, "num_min_like": int(min_like.size), "num_max_like": int(max_like.size),
        "min_candidates_scored": len(min_rows), "max_candidates_scored": len(max_rows),
        "selected_min": sum(r["polarity"] == "min" for r in selected),
        "selected_max": sum(r["polarity"] == "max" for r in selected),
        "max_min_height": max_min_height, "max_max_height": max_max_height,
        "avg_min_like_area": avg_min_area, "avg_max_like_area": avg_max_area,
        "dark_alt_max": dark_alt_max, "bright_alt_min": bright_alt_min,
        "scoring_seconds": time.time() - t1,
    }
    log(f"[tos] scored MIN={len(min_rows)} MAX={len(max_rows)} -> selected "
        f"{info['selected_min']} MIN + {info['selected_max']} MAX in {info['scoring_seconds']:.1f}s")
    return tree, selected, info


# ============================================================ spatial diagnostics (reconstructed e10, BG9/BG12)
def exclusive_support_voxels(tree, selected_nodes: np.ndarray) -> np.ndarray:
    """Return, for each leaf voxel, the id of its nearest selected ancestor (or -1)."""
    n = tree.num_vertices()
    selected_mask = np.zeros(n, dtype=bool)
    selected_mask[selected_nodes] = True
    owner = hg.propagate_sequential(tree, np.arange(n, dtype=np.int64), ~selected_mask)
    owner = np.asarray(owner)[: tree.num_leaves()].copy()
    owner[~selected_mask[owner]] = -1
    return owner


def spatialize(tree, candidates: list[dict], shape, internal: np.ndarray, cfg: dict[str, Any]) -> pd.DataFrame:
    labels_internal, _ = ndi.label(internal, structure=CONNECTIVITY_6)
    n_leaves = tree.num_leaves()
    support_nodes = np.asarray([c["node"] for c in candidates
                                if cfg["support_set"] == "all" or c["polarity"] == "min"], dtype=np.int64)
    owner = exclusive_support_voxels(tree, support_nodes)
    counts = np.bincount(owner[owner >= 0], minlength=tree.num_vertices())
    limits = np.asarray(shape) - 1
    rows = []
    for cand in candidates:
        node = cand["node"]
        _, node_map = tree.sub_tree(node)
        leaves = np.asarray(node_map)
        leaves = leaves[leaves < n_leaves]
        row = dict(cand)
        if leaves.size == 0:
            row.update(has_reconstructed_region=False, point_z=np.nan, point_y=np.nan, point_x=np.nan,
                       region_voxels=0, region_full_voxels=0, region_touches_volume_face=False,
                       region_internal_overlap_fraction=0.0, region_dominant_internal_purity=0.0,
                       region_internal_components=0)
            rows.append(row)
            continue
        coords = np.stack(np.unravel_index(leaves, shape), axis=1)
        lo, hi = coords.min(axis=0), coords.max(axis=0)
        if cfg["point_mode"] == "medial":
            box = np.zeros(tuple(hi - lo + 3), dtype=bool)
            box[tuple((coords - lo + 1).T)] = True
            edt = ndi.distance_transform_edt(box)
            local = np.unravel_index(int(np.argmax(edt)), box.shape)
            point = np.asarray(local) + lo - 1
        elif cfg["point_mode"] == "centroid_snapped":
            centroid = coords.mean(axis=0)
            point = coords[int(np.argmin(((coords - centroid) ** 2).sum(axis=1)))]
        else:
            raise ValueError(f"Unknown tos.point_mode {cfg['point_mode']!r}")
        inside = internal[tuple(coords.T)]
        comps = labels_internal[tuple(coords.T)]
        comps = comps[comps > 0]
        if comps.size:
            _, comp_counts = np.unique(comps, return_counts=True)
            purity = float(comp_counts.max()) / len(coords)
            n_comps = int(comp_counts.size)
        else:
            purity, n_comps = 0.0, 0
        row.update(
            has_reconstructed_region=True,
            point_z=int(point[0]), point_y=int(point[1]), point_x=int(point[2]),
            region_voxels=int(counts[node]),
            region_full_voxels=int(len(coords)),
            region_touches_volume_face=bool((lo <= 0).any() or (hi >= limits).any()),
            region_internal_overlap_fraction=float(inside.mean()),
            region_dominant_internal_purity=purity,
            region_internal_components=n_comps,
        )
        rows.append(row)
    return pd.DataFrame(rows)


# ============================================================ markers (exact port of e11_build_markers_ablation)
def normalize_score(series: pd.Series) -> np.ndarray:
    numeric = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float).copy()
    finite = np.isfinite(numeric)
    if not finite.any():
        return np.full(len(numeric), 0.5, dtype=float)
    numeric[~finite] = float(np.nanmedian(numeric[finite]))
    low, high = float(np.min(numeric)), float(np.max(numeric))
    if high <= low:
        return np.full(len(numeric), 0.5, dtype=float)
    return (numeric - low) / (high - low)


def greedy_nms(frame: pd.DataFrame, radius: float) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    if radius <= 0:
        raise ValueError("nms-radius must be positive.")
    physical = frame[["projected_z", "projected_y", "projected_x"]].to_numpy(dtype=np.int64).astype(float)
    priority = frame["priority"].to_numpy(dtype=float)
    order = np.lexsort((frame["source_row_index"].to_numpy(dtype=np.int64), -priority))
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


def build_markers(frame: pd.DataFrame, internal: np.ndarray, cfg: dict[str, Any]):
    labels, count = ndi.label(internal, structure=CONNECTIVITY_6)
    if count == 0:
        raise ValueError("internal_dark mask is empty.")
    edt = ndi.distance_transform_edt(internal)
    positions = ndi.maximum_position(edt, labels=labels, index=np.arange(1, count + 1, dtype=np.int32))
    max_edt = np.asarray([edt[tuple(p)] for p in np.asarray(positions, dtype=np.int32)], dtype=np.float32)

    coords_float = np.column_stack([pd.to_numeric(frame[c], errors="coerce") for c in ("point_z", "point_y", "point_x")])
    has_region = frame["has_reconstructed_region"].astype(bool).to_numpy()
    valid = has_region & np.isfinite(coords_float).all(axis=1)
    frame = frame.loc[valid].copy().reset_index(drop=True)
    coords = np.rint(coords_float[valid]).astype(np.int64)
    inside = ((coords >= 0) & (coords < np.asarray(internal.shape))).all(axis=1)
    frame = frame.loc[inside].copy().reset_index(drop=True)
    coords = coords[inside]

    polarity = frame["polarity"].astype(str).str.lower().str.strip()
    base = pd.DataFrame({
        "source_row_index": np.arange(len(frame), dtype=np.int64),
        "node": pd.to_numeric(frame["node"], errors="raise").astype(np.int64),
        "polarity": polarity.to_numpy(),
        "raw_score": pd.to_numeric(frame["score"], errors="coerce"),
        "score_norm": normalize_score(frame["score"]),
        "original_z": coords[:, 0], "original_y": coords[:, 1], "original_x": coords[:, 2],
    })
    base["point_internal_component"] = labels[tuple(coords[:, a] for a in range(3))].astype(np.int64)
    if cfg["polarity"] == "dual":
        polarity_mask = base["polarity"].isin(["min", "max"])
    else:
        polarity_mask = base["polarity"].eq(cfg["polarity"])
    selected = base.loc[polarity_mask & (base["point_internal_component"] > 0)].copy()
    if selected.empty:
        raise RuntimeError(f"No {cfg['polarity'].upper()} candidate falls inside internal_dark.")
    component_ids = selected["point_internal_component"].to_numpy(dtype=np.int64)
    edt_value = max_edt[component_ids - 1]
    maximum = float(max_edt.max())
    edt_norm = np.log1p(edt_value) / math.log1p(maximum) if maximum > 0 else np.zeros(len(selected))
    min_bonus = selected["polarity"].eq("min").to_numpy(dtype=float) * float(cfg["min_priority_bonus"])
    selected["mapped_internal_component"] = component_ids
    selected["component_max_edt"] = edt_value
    selected["projection_distance"] = 0.0
    selected["priority"] = (selected["score_norm"].to_numpy(dtype=float)
                            + float(cfg["edt_priority_weight"]) * edt_norm + min_bonus)
    selected = (selected.sort_values(["mapped_internal_component", "priority", "score_norm", "source_row_index"],
                                     ascending=[True, False, False, True])
                .drop_duplicates("mapped_internal_component", keep="first").reset_index(drop=True))
    selected["projected_z"] = selected["original_z"].astype(np.int64)
    selected["projected_y"] = selected["original_y"].astype(np.int64)
    selected["projected_x"] = selected["original_x"].astype(np.int64)
    before_nms = int(len(selected))
    markers = greedy_nms(selected, float(cfg["nms_radius"]))
    volume = np.zeros(internal.shape, dtype=np.uint16 if len(markers) <= 65535 else np.uint32)
    for marker_id, point in enumerate(markers[["projected_z", "projected_y", "projected_x"]].to_numpy(np.int64), 1):
        volume[tuple(point)] = marker_id
    info = {"input_spatial_candidates": int(len(base)),
            "candidates_inside_internal": int((polarity_mask & (base["point_internal_component"] > 0)).sum()),
            "markers_before_nms": before_nms, "marker_count": int(len(markers)),
            "internal_components": int(count)}
    return markers, volume, info


# ============================================================ stage runner
def run(selector_u8: np.ndarray, config: dict[str, Any], out_dir: Path, prefix: str, log=print) -> dict[str, Any]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    selector_u8 = robust_u8(selector_u8)
    internal, boundary_info = internal_dark_mask(selector_u8, config["boundary"])
    internal_path = write3d(out_dir / f"{prefix}_internal_dark_candidates.tif", internal.astype(np.uint8) * 255)
    log(f"[boundary] dark_thr={boundary_info['dark_threshold']} internal voxels={boundary_info['internal_dark_voxels']}")

    tree, candidates, tos_info = score_candidates(selector_u8, config["tos"], log=log)
    cand_path = out_dir / f"{prefix}_tosconop_autodual_candidates.csv"
    pd.DataFrame(candidates).assign(rank=lambda f: np.arange(1, len(f) + 1)).to_csv(cand_path, index=False)

    t0 = time.time()
    diagnostics = spatialize(tree, candidates, selector_u8.shape, internal, config["tos"])
    if config["boundary"].get("reject_node_touching_face", False):
        diagnostics.loc[diagnostics["region_touches_volume_face"], "has_reconstructed_region"] = False
    diag_path = out_dir / "candidate_diagnostics_no_gt.csv"
    diagnostics.to_csv(diag_path, index=False)
    log(f"[spatial] {len(diagnostics)} candidates spatialized in {time.time() - t0:.1f}s")

    # Port (cross-check; always computed)
    markers, volume, marker_info = build_markers(diagnostics, internal, config["markers"])
    markers.to_csv(out_dir / "markers_port.csv", index=False)
    port_path = write3d(out_dir / "markers_port_uint16.tif", volume)
    markers_csv, marker_path = out_dir / "markers_port.csv", port_path
    if config["markers"].get("implementation", "port") == "e11_verbatim":
        from .recovered import e11_markers
        verbatim = e11_markers(config, diag_path, internal_path, out_dir / "markers_e11")
        same = bool(np.array_equal(read3d(verbatim["markers_tif"]), volume))
        marker_info = dict(marker_info, **verbatim, port_identical_to_verbatim=same)
        if not same:
            log("[markers] WARNING: port differs from verbatim E11 builder; using verbatim")
        markers_csv, marker_path = Path(verbatim["markers_csv"]), Path(verbatim["markers_tif"])
    log(f"[markers] {marker_info['markers_before_nms']} consolidated -> {marker_info['marker_count']} after NMS "
        f"({config['markers'].get('implementation', 'port')})")
    info = {"boundary": boundary_info, "tos": tos_info, "markers": marker_info,
            "internal_dark": str(internal_path), "candidates_csv": str(cand_path),
            "diagnostics_csv": str(diag_path), "markers_csv": str(markers_csv),
            "markers_tif": str(marker_path)}
    write_json(out_dir / "tos_select.json", info)
    return info


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="batos-paper-v2 stage 2: ToS selection and markers")
    parser.add_argument("selector", type=Path, help="*_flatten_lambda1p0_8u.tif")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--config", type=Path, default=Path("configs/batos_paper_v2.yaml"))
    args = parser.parse_args(argv)
    run(read3d(args.selector, "selector"), load_config(args.config), args.out_dir, args.prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
