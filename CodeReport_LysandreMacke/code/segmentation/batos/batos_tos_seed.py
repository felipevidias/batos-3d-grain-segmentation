#!/usr/bin/env python3
"""
Autodual ToSConOp seed extraction for the Lysandre/Macke watershed pipeline.

This script keeps the successful Macke downstream protocol:
    binary core mask -> 3D connected components -> centroid markers -> watershed

but replaces only the Min-tree core extraction by a self-dual Tree of Shapes strategy:
    MIN-like dark nodes  = candidate grain cores
    MAX-like ancestors   = autodual support/ring evidence

The output seed is intentionally a binary core mask compatible with watershed.py.
"""

import argparse
import csv
from pathlib import Path
import warnings

import numpy as np
import tifffile
from scipy import ndimage
from skimage.measure import label, regionprops
import skimage.morphology as morph

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import higra as hg


def read_tiff_stack(path: str) -> np.ndarray:
    arr = tifffile.imread(path)
    arr = np.asarray(arr)
    if arr.ndim == 2:
        arr = arr[np.newaxis, :, :]
    if arr.ndim == 4 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    if arr.ndim != 3:
        raise ValueError(f"Expected 3D TIFF stack, got shape {arr.shape}")
    return arr


def robust_to_uint8(volume: np.ndarray, p_low: float = 0.5, p_high: float = 99.5) -> np.ndarray:
    if volume.dtype == np.uint8:
        return volume.copy()
    vol = volume.astype(np.float32, copy=False)
    finite = np.isfinite(vol)
    nz = vol[(vol > 0) & finite]
    if nz.size == 0:
        nz = vol[finite]
    if nz.size == 0:
        return np.zeros_like(vol, dtype=np.uint8)
    lo = np.percentile(nz, p_low)
    hi = np.percentile(nz, p_high)
    if hi <= lo:
        hi = lo + 1.0
    out = (vol - lo) * (255.0 / (hi - lo))
    return np.clip(out, 0, 255).astype(np.uint8)


def flatten_volume_u8(volume_u8: np.ndarray,
                      sigma_xy: float = 35.0,
                      sigma_z: float = 0.0,
                      p_low: float = 0.5,
                      p_high: float = 99.5) -> np.ndarray:
    vol = volume_u8.astype(np.float32)
    if sigma_z > 0:
        background = ndimage.gaussian_filter(vol, sigma=(sigma_z, sigma_xy, sigma_xy))
    else:
        background = np.empty_like(vol)
        for z in range(vol.shape[0]):
            background[z] = ndimage.gaussian_filter(vol[z], sigma=sigma_xy)
    nonzero_bg = background[background > 0]
    bg_ref = float(np.median(nonzero_bg)) if nonzero_bg.size else float(np.median(background))
    flat = vol - background + bg_ref
    flat = np.clip(flat, 0, None)
    nz = flat[flat > 0]
    if nz.size == 0:
        return np.zeros_like(volume_u8, dtype=np.uint8)
    lo = np.percentile(nz, p_low)
    hi = np.percentile(nz, p_high)
    if hi <= lo:
        hi = lo + 1.0
    flat8 = ((flat - lo) * (255.0 / (hi - lo)))
    return np.clip(flat8, 0, 255).astype(np.uint8)


def clamp01(x: float) -> float:
    return float(max(0.0, min(1.0, x)))


def node_is_min_like(node: int, parents: np.ndarray, altitudes: np.ndarray, root: int) -> bool:
    if node == root:
        return False
    parent = int(parents[node])
    if parent == node or parent < 0:
        return False
    return int(altitudes[node]) < int(altitudes[parent])


def node_is_max_like(node: int, parents: np.ndarray, altitudes: np.ndarray, root: int) -> bool:
    if node == root:
        return False
    parent = int(parents[node])
    if parent == node or parent < 0:
        return False
    return int(altitudes[node]) > int(altitudes[parent])


def nearest_max_support(node: int,
                        parents: np.ndarray,
                        altitudes: np.ndarray,
                        area: np.ndarray,
                        root: int,
                        support_area_min: int,
                        support_area_max: int,
                        max_steps: int):
    cur = int(parents[node])
    steps = 1
    best = -1
    while cur != root and cur >= 0 and steps <= max_steps:
        if node_is_max_like(cur, parents, altitudes, root):
            if support_area_min <= int(area[cur]) <= support_area_max:
                best = cur
                break
        nxt = int(parents[cur])
        if nxt == cur:
            break
        cur = nxt
        steps += 1
    return best, steps


def build_autodual_candidates(image_u8: np.ndarray, args):
    print("[ToSConOp autodual] Building 3D Tree of Shapes with Higra...")
    tree, altitudes = hg.component_tree_tree_of_shapes_image3d(image_u8)
    altitudes = altitudes.astype(np.int64)
    parents = tree.parents().astype(np.int64)
    root = int(tree.root())

    print("[ToSConOp autodual] Computing attributes...")
    area = hg.attribute_area(tree).astype(np.int64)
    height = hg.attribute_height(tree, altitudes).astype(np.float64)

    min_like = []
    max_like = []
    for node in range(len(parents)):
        if node_is_min_like(node, parents, altitudes, root):
            min_like.append(node)
        elif node_is_max_like(node, parents, altitudes, root):
            max_like.append(node)

    min_like = np.asarray(min_like, dtype=np.int64)
    max_like = np.asarray(max_like, dtype=np.int64)

    max_min_height = float(np.max(height[min_like])) if min_like.size else 1.0
    max_max_height = float(np.max(height[max_like])) if max_like.size else 1.0
    avg_min_area = float(np.average(area[min_like])) if min_like.size else 1.0
    nz = image_u8[image_u8 > 0]
    dark_alt_max = float(np.percentile(nz, args.dark_percentile_max)) if nz.size else 255.0

    candidates = []
    # Broad prefilter first; scoring is where autodual evidence enters.
    for node in min_like:
        node = int(node)
        parent = int(parents[node])
        if parent < 0 or parent == node:
            continue
        node_area = int(area[node])
        if node_area < args.area_min or node_area > args.area_max:
            continue
        if node_area > args.area_factor * avg_min_area:
            continue
        if int(altitudes[node]) > dark_alt_max:
            continue
        delta_parent = int(altitudes[parent]) - int(altitudes[node])
        if delta_parent < args.delta_parent_min:
            continue

        support, support_steps = nearest_max_support(
            node, parents, altitudes, area, root,
            support_area_min=args.support_area_min,
            support_area_max=args.support_area_max,
            max_steps=args.support_max_steps,
        )
        has_support = support >= 0
        if args.require_max_support and not has_support:
            continue

        h_score = clamp01(float(height[node]) / max(1e-9, args.height_factor * max_min_height))
        delta_score = clamp01(delta_parent / max(1.0, args.delta_norm))
        # Macke-like preference: small but not degenerate core.
        log_area = np.log(float(node_area) + 1.0)
        log_min = np.log(float(max(1, args.area_min)) + 1.0)
        log_max = np.log(float(max(args.area_min + 1, args.area_max)) + 1.0)
        area_pos = (log_area - log_min) / max(1e-9, log_max - log_min)
        area_score = clamp01(1.0 - abs(area_pos - args.area_preferred_position) / max(1e-9, args.area_preferred_width))
        dark_score = clamp01((dark_alt_max - float(altitudes[node])) / max(1.0, dark_alt_max))

        support_score = 0.0
        support_area = 0
        support_height = 0.0
        support_contrast = 0.0
        if has_support:
            support_area = int(area[support])
            support_height = float(height[support])
            support_contrast = abs(float(altitudes[support]) - float(altitudes[node]))
            support_height_score = clamp01(support_height / max(1e-9, 0.10 * max_max_height))
            ratio = float(support_area) / max(1.0, float(node_area))
            ratio_score = clamp01(1.0 - abs(np.log(ratio + 1e-9) - np.log(args.support_area_ratio_preferred)) /
                                  max(1e-9, np.log(args.support_area_ratio_width)))
            contrast_score = clamp01(support_contrast / max(1.0, args.support_contrast_norm))
            depth_score = clamp01(1.0 - (support_steps - 1) / max(1.0, args.support_max_steps))
            support_score = 0.35 * support_height_score + 0.30 * ratio_score + 0.25 * contrast_score + 0.10 * depth_score

        score = (
            args.w_height * h_score +
            args.w_delta * delta_score +
            args.w_area * area_score +
            args.w_dark * dark_score +
            args.w_support * support_score
        )

        if score < args.min_score:
            continue

        candidates.append({
            "node": node,
            "score": score,
            "area": node_area,
            "height": float(height[node]),
            "alt": int(altitudes[node]),
            "parent_alt": int(altitudes[parent]),
            "delta_parent": delta_parent,
            "support": int(support),
            "support_area": int(support_area),
            "support_height": float(support_height),
            "support_contrast": float(support_contrast),
            "support_steps": int(support_steps),
            "h_score": h_score,
            "delta_score": delta_score,
            "area_score": area_score,
            "dark_score": dark_score,
            "support_score": support_score,
        })

    candidates.sort(key=lambda c: (c["score"], c["height"], -c["area"]), reverse=True)

    if args.max_selected_nodes > 0:
        candidates = candidates[:args.max_selected_nodes]

    info = {
        "num_tree_nodes": int(len(parents)),
        "num_min_like_nodes": int(len(min_like)),
        "num_max_like_nodes": int(len(max_like)),
        "num_candidates_after_autodual_scoring": int(len(candidates)),
        "max_min_height": float(max_min_height),
        "max_max_height": float(max_max_height),
        "avg_min_like_area": float(avg_min_area),
        "dark_alt_max": float(dark_alt_max),
    }
    print(f"[ToSConOp autodual] nodes={len(parents)} min_like={len(min_like)} max_like={len(max_like)} candidates={len(candidates)}")
    print(f"[ToSConOp autodual] max_min_height={max_min_height:.3f} avg_min_area={avg_min_area:.3f} dark_alt_max={dark_alt_max:.3f}")
    return tree, candidates, info


def reconstruct_seed_from_candidates(tree, candidates, shape):
    parents_len = len(tree.parents())
    node_values = np.zeros(parents_len, dtype=np.uint8)
    for c in candidates:
        node = int(c["node"])
        try:
            _, subnodes = tree.sub_tree(node)
            node_values[subnodes] = 1
        except Exception:
            node_values[node] = 1
    seed = hg.reconstruct_leaf_data(tree, node_values)
    seed = np.asarray(seed).reshape(shape)
    return (seed > 0).astype(np.uint8) * 255


def cleanup_seed(seed: np.ndarray, args):
    binary = seed > 0
    if args.cleanup_open_radius > 0:
        binary = morph.opening(binary, morph.ball(args.cleanup_open_radius))
    if args.cleanup_close_radius > 0:
        binary = morph.closing(binary, morph.ball(args.cleanup_close_radius))

    lab = label(binary, background=0, connectivity=1)
    out = np.zeros_like(seed, dtype=np.uint8)
    centroids = np.zeros_like(seed, dtype=np.uint16)

    kept = 0
    sizes = []
    for region in regionprops(lab):
        n = int(region.area)
        if n < args.component_min or n > args.component_max:
            continue
        kept += 1
        sizes.append(n)
        coords = region.coords
        out[coords[:, 0], coords[:, 1], coords[:, 2]] = 255
        cz, cy, cx = np.round(region.centroid).astype(int)
        cz = int(np.clip(cz, 0, seed.shape[0] - 1))
        cy = int(np.clip(cy, 0, seed.shape[1] - 1))
        cx = int(np.clip(cx, 0, seed.shape[2] - 1))
        if kept <= np.iinfo(np.uint16).max:
            centroids[cz, cy, cx] = kept

    info = {
        "num_components_kept": int(kept),
        "seed_voxels": int(np.count_nonzero(out)),
        "seed_mean_0_255": float(np.mean(out)),
        "component_size_min": int(min(sizes)) if sizes else 0,
        "component_size_max": int(max(sizes)) if sizes else 0,
        "component_size_mean": float(np.mean(sizes)) if sizes else 0.0,
    }
    return out, centroids, info


def write_report(path: Path, args, tree_info, comp_info):
    with open(path, "w", encoding="utf-8") as f:
        f.write("Autodual ToSConOp seed extraction report\n")
        f.write("========================================\n\n")
        f.write("Concept: MIN-like dark ToS cores validated/scored by MAX-like ancestor support.\n")
        f.write("Downstream protocol remains Macke-style centroid watershed.\n\n")
        f.write("Command parameters:\n")
        for k, v in vars(args).items():
            f.write(f"  {k}: {v}\n")
        f.write("\nTree statistics:\n")
        for k, v in tree_info.items():
            f.write(f"  {k}: {v}\n")
        f.write("\nSeed/component statistics:\n")
        for k, v in comp_info.items():
            f.write(f"  {k}: {v}\n")


def write_candidate_csv(path: Path, candidates):
    fields = [
        "rank", "node", "score", "area", "height", "alt", "parent_alt", "delta_parent",
        "support", "support_area", "support_height", "support_contrast", "support_steps",
        "h_score", "delta_score", "area_score", "dark_score", "support_score",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for i, c in enumerate(candidates, start=1):
            row = {k: c.get(k, "") for k in fields}
            row["rank"] = i
            writer.writerow(row)


def main():
    ap = argparse.ArgumentParser(description="Autodual ToSConOp seed extraction for Macke-style watershed.")
    ap.add_argument("image")
    ap.add_argument("--out-dir", default="results_tosconop_autodual_seed")
    ap.add_argument("--flatten", action="store_true")
    ap.add_argument("--flatten-sigma-xy", type=float, default=35.0)
    ap.add_argument("--flatten-sigma-z", type=float, default=0.0)

    # Broad Macke-like core constraints.
    ap.add_argument("--height-factor", type=float, default=0.003)
    ap.add_argument("--area-factor", type=float, default=80.0)
    ap.add_argument("--area-min", type=int, default=2)
    ap.add_argument("--area-max", type=int, default=16000)
    ap.add_argument("--delta-parent-min", type=int, default=0)
    ap.add_argument("--dark-percentile-max", type=float, default=99.0)

    # Autodual support constraints.
    ap.add_argument("--support-area-min", type=int, default=30)
    ap.add_argument("--support-area-max", type=int, default=90000)
    ap.add_argument("--support-max-steps", type=int, default=28)
    ap.add_argument("--require-max-support", action="store_true")
    ap.add_argument("--support-area-ratio-preferred", type=float, default=18.0)
    ap.add_argument("--support-area-ratio-width", type=float, default=80.0)
    ap.add_argument("--support-contrast-norm", type=float, default=80.0)

    # Candidate score.
    ap.add_argument("--delta-norm", type=float, default=20.0)
    ap.add_argument("--area-preferred-position", type=float, default=0.42)
    ap.add_argument("--area-preferred-width", type=float, default=0.55)
    ap.add_argument("--w-height", type=float, default=0.22)
    ap.add_argument("--w-delta", type=float, default=0.12)
    ap.add_argument("--w-area", type=float, default=0.18)
    ap.add_argument("--w-dark", type=float, default=0.18)
    ap.add_argument("--w-support", type=float, default=0.30)
    ap.add_argument("--min-score", type=float, default=0.22)
    ap.add_argument("--max-selected-nodes", type=int, default=3500)

    # Cleanup / components.
    ap.add_argument("--cleanup-open-radius", type=int, default=0)
    ap.add_argument("--cleanup-close-radius", type=int, default=0)
    ap.add_argument("--component-min", type=int, default=2)
    ap.add_argument("--component-max", type=int, default=12000)
    args = ap.parse_args()

    image_path = Path(args.image)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    base = image_path.stem

    print(f"[ToSConOp autodual] Reading {image_path}")
    raw = read_tiff_stack(str(image_path))
    image_u8 = robust_to_uint8(raw)
    if args.flatten:
        print("[ToSConOp autodual] Applying flattening before ToS...")
        tos_input = flatten_volume_u8(image_u8, sigma_xy=args.flatten_sigma_xy, sigma_z=args.flatten_sigma_z)
        tifffile.imwrite(out_dir / f"{base}_tosconop_autodual_flattened_8u.tif", tos_input, imagej=True)
    else:
        tos_input = image_u8
        tifffile.imwrite(out_dir / f"{base}_tosconop_autodual_input_8u.tif", tos_input, imagej=True)

    tree, candidates, tree_info = build_autodual_candidates(tos_input, args)
    write_candidate_csv(out_dir / f"{base}_tosconop_autodual_candidates.csv", candidates)

    seed = reconstruct_seed_from_candidates(tree, candidates, tos_input.shape)
    seed_clean, centroids, comp_info = cleanup_seed(seed, args)

    seed_path = out_dir / f"{base}_tosconop_autodual_seed_raw.tif"
    centroid_path = out_dir / f"{base}_tosconop_autodual_seed_centroids.tif"
    report_path = out_dir / f"{base}_tosconop_autodual_seed_report.txt"
    tifffile.imwrite(seed_path, seed_clean.astype(np.uint8), imagej=True)
    tifffile.imwrite(centroid_path, centroids.astype(np.uint16), imagej=True)
    write_report(report_path, args, tree_info, comp_info)

    print(f"[ToSConOp autodual] Wrote seed mask    : {seed_path}")
    print(f"[ToSConOp autodual] Wrote centroids    : {centroid_path}")
    print(f"[ToSConOp autodual] Wrote candidates   : {out_dir / (base + '_tosconop_autodual_candidates.csv')}")
    print(f"[ToSConOp autodual] Wrote report       : {report_path}")
    print(f"[ToSConOp autodual] selected nodes     = {len(candidates)}")
    print(f"[ToSConOp autodual] seed voxels        = {comp_info['seed_voxels']}")
    print(f"[ToSConOp autodual] core components    = {comp_info['num_components_kept']}")
    print(f"[ToSConOp autodual] seed mean 0..255   = {comp_info['seed_mean_0_255']:.6f}")


if __name__ == "__main__":
    main()
