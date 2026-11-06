"""Phase 5c: extinction-based ToS marker selection (candidate replacement for the frozen weighted score).

Drop-in replacement for the frozen true_dual selector (batos_tos_seed_true_dual.py, K=3600): it writes a candidate
CSV with the same columns and file name, so the frozen e10 spatialisation, e11 marker builder, EDT-Otsu control,
E23c gate, fixed-mask watershed and evaluate_v2 run unchanged downstream. The frozen code is not modified.

Algorithm (one crop, flattened selector image, uint8):
 1. Tree of Shapes with Higra, exactly as the frozen selector and e10 build it
    (hg.component_tree_tree_of_shapes_image3d on the uint8 image; e10 re-checks node ids/altitudes).
 2. Dark extrema M = Higra extrema (hg.attribute_extrema) that are MIN-like (alt(n) < alt(parent(n))).
 3. Boundary rule: a node "touches a face" if its shape contains a voxel on any of the 6 crop faces. Extrema that
    touch a face are rejected (pool = remaining dark extrema).
 4. Extinction values on the ToS. Higra's attribute_extinction_value/attribute_dynamics assume monotone altitudes
    (min/max trees) and do not apply to a ToS, so the classic branch rule is computed directly for an attribute A
    that is increasing on the tree: at every node, among the children whose subtree contains a dark extremum, the
    child with the largest A wins (ties: smallest node index); every other such child is the top of a branch that
    dies there. The extinction value of a dark extremum m is A(top of m's branch) (A(root) for the global winner).
      volume   : A = hg.attribute_volume (sum over the subtree of area * |alt - alt(parent)|)
      dynamics : A = max - min grey level inside the shape (ToS analogue of the height/dynamics extinction)
      product  : ext_volume(m) * ext_dynamics(m) (each from its own competition)
 5. Rank the pool by extinction value (ties: darker extremum first, then node index); keep the top K = 3600.
 6. Representative node of each selected extremum m: the highest node on m's own branch (m up to its branch top)
    whose shape contains no other selected extremum and touches no face, restricted to MIN-like nodes. These
    shapes are pairwise disjoint (each contains exactly one selected extremum), so no hierarchical suppression is
    needed.
 7. CSV `score` = rank score (K - rank + 1) / K in (0, 1]: the frozen e11 builder min-max normalises `score` and
    mixes it with an EDT term; a monotone, parameter-free rank keeps that mix meaningful for heavy-tailed
    extinction values. The raw value is written to `extinction`.

Free parameters of this stage: K (3600, as frozen) and the variant (categorical).
Frozen weighted selector: 17 swept parameters (configs/sensitivity_global.yaml, group `selector`).
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import higra as hg
import numpy as np
import tifffile

VARIANTS = ("volume", "dynamics", "product")
FIELDS = ["rank", "node", "polarity", "score", "area", "height", "polarity_height", "alt", "parent_alt",
          "delta_parent", "support", "support_area", "support_height", "support_contrast", "support_steps",
          "h_score", "delta_score", "area_score", "dark_score", "intensity_score", "support_score",
          "extinction", "extinction_volume", "extinction_dynamics", "extremum_node", "extremum_alt"]


def branch_extinction(parents: np.ndarray, root: int, attribute: np.ndarray, has_dark: np.ndarray, tree):
    """Per-node extinction value of the branch the node belongs to, and the 'winner child' mask."""
    n = len(parents)
    child = np.arange(n, dtype=np.int64)
    child = child[child != root]
    par = parents[child]
    key = np.where(has_dark[child], attribute[child].astype(np.float64), -1.0)
    best = np.full(n, -np.inf)
    np.maximum.at(best, par, key)
    win = has_dark[child] & (key == best[par])
    first = np.full(n, n, dtype=np.int64)
    np.minimum.at(first, par[win], child[win])
    win &= child == first[par]
    is_winner = np.zeros(n, dtype=bool)
    is_winner[child[win]] = True
    ext = hg.propagate_sequential(tree, attribute.astype(np.float64), is_winner)
    return ext, is_winner


def select(image: np.ndarray, variant: str, k: int) -> tuple[list[dict], dict]:
    if variant not in VARIANTS:
        raise ValueError(variant)
    t0 = time.time()
    tree, altitudes = hg.component_tree_tree_of_shapes_image3d(image)
    alt = np.asarray(altitudes, dtype=np.int64)
    parents = tree.parents().astype(np.int64)
    root, nl, n = int(tree.root()), int(tree.num_leaves()), len(parents)
    if nl != image.size:
        raise RuntimeError(f"leaves={nl}, voxels={image.size}")
    area = hg.attribute_area(tree).astype(np.int64)
    min_like = alt < alt[parents]
    min_like[root] = False
    dark = np.asarray(hg.attribute_extrema(tree, alt), dtype=bool) & min_like

    face_leaf = np.zeros(image.shape, dtype=np.uint8)
    face_leaf[0], face_leaf[-1] = 1, 1
    face_leaf[:, 0], face_leaf[:, -1] = 1, 1
    face_leaf[:, :, 0], face_leaf[:, :, -1] = 1, 1
    touches = hg.accumulate_sequential(tree, face_leaf.ravel(), hg.Accumulators.max).astype(bool)
    has_dark = hg.accumulate_and_max_sequential(tree, dark.astype(np.uint8), np.zeros(nl, np.uint8),
                                                hg.Accumulators.max).astype(bool)

    ext = {}
    winners = {}
    if variant in ("volume", "product"):
        vol = hg.attribute_volume(tree, alt.astype(np.float64), area.astype(np.float64))
        ext["volume"], winners["volume"] = branch_extinction(parents, root, vol, has_dark, tree)
        del vol
    if variant in ("dynamics", "product"):
        leaf_alt = alt[:nl]
        rng = (hg.accumulate_sequential(tree, leaf_alt, hg.Accumulators.max)
               - hg.accumulate_sequential(tree, leaf_alt, hg.Accumulators.min)).astype(np.float64)
        ext["dynamics"], winners["dynamics"] = branch_extinction(parents, root, rng, has_dark, tree)
        del rng
    if variant == "product":
        value = ext["volume"] * ext["dynamics"]
        # branch membership for the representative walk: a node stays on m's branch only if it wins in BOTH
        # competitions (conservative: the region of m never extends past either extinction point)
        is_winner = winners["volume"] & winners["dynamics"]
    else:
        value = ext[variant]
        is_winner = winners[variant]

    dark_nodes = np.flatnonzero(dark)
    pool = dark_nodes[~touches[dark_nodes]]
    order = np.lexsort((pool, alt[pool], -value[pool]))
    chosen = pool[order[:k]]

    selected = np.zeros(n, dtype=np.uint8)
    selected[chosen] = 1
    count = hg.accumulate_and_add_sequential(tree, selected.astype(np.int64), np.zeros(nl, np.int64),
                                             hg.Accumulators.sum)
    leaf_rng = alt[:nl]
    hmax = hg.accumulate_sequential(tree, leaf_rng, hg.Accumulators.max)
    hmin = hg.accumulate_sequential(tree, leaf_rng, hg.Accumulators.min)

    rows, depth_steps = [], []
    for rank, m in enumerate(chosen.tolist(), start=1):
        cur, rep, steps = m, m, 0
        while cur != root and is_winner[cur]:
            p = int(parents[cur])
            if touches[p] or count[p] > 1:
                break
            cur = p
            steps += 1
            if min_like[cur]:
                rep = cur
        depth_steps.append(steps)
        pa = int(parents[rep])
        h = float(hmax[rep] - hmin[rep])
        rows.append({
            "rank": rank, "node": rep, "polarity": "min", "score": (k - rank + 1) / k if k > 0 else 1.0,
            "area": int(area[rep]), "height": h, "polarity_height": h, "alt": int(alt[rep]),
            "parent_alt": int(alt[pa]), "delta_parent": int(alt[pa] - alt[rep]), "support": -1, "support_area": 0,
            "support_height": 0.0, "support_contrast": 0.0, "support_steps": 0, "h_score": "", "delta_score": "",
            "area_score": "", "dark_score": "", "intensity_score": "", "support_score": 0.0,
            "extinction": float(value[m]),
            "extinction_volume": float(ext["volume"][m]) if "volume" in ext else "",
            "extinction_dynamics": float(ext["dynamics"][m]) if "dynamics" in ext else "",
            "extremum_node": int(m), "extremum_alt": int(alt[m])})
    nodes = [r["node"] for r in rows]
    if len(set(nodes)) != len(nodes):
        raise RuntimeError("representative nodes are not unique")
    areas = np.array([r["area"] for r in rows]) if rows else np.zeros(0)
    info = {"variant": variant, "k": k, "num_tree_nodes": n, "num_leaves": nl,
            "num_min_like_nodes": int(min_like.sum()), "num_dark_extrema": int(dark_nodes.size),
            "num_dark_extrema_touching_face": int(dark_nodes.size - pool.size), "pool_size": int(pool.size),
            "num_selected": len(rows), "representative_area_median": float(np.median(areas)) if rows else None,
            "representative_area_q75": float(np.percentile(areas, 75)) if rows else None,
            "representative_area_max": int(areas.max()) if rows else None,
            "representative_steps_median": float(np.median(depth_steps)) if rows else None,
            "extinction_kth": float(value[chosen[-1]]) if rows else None, "seconds": round(time.time() - t0, 1)}
    return rows, info


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Phase 5c extinction-based ToS selector (true_dual replacement)")
    ap.add_argument("image", type=Path)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--variant", choices=VARIANTS, required=True)
    ap.add_argument("--k", type=int, default=3600)
    args = ap.parse_args(argv)
    image = np.squeeze(np.asarray(tifffile.imread(args.image)))
    if image.ndim != 3 or image.dtype != np.uint8:
        raise RuntimeError(f"selector image must be 3D uint8, got {image.shape} {image.dtype}")
    rows, info = select(image, args.variant, args.k)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    base = args.image.stem
    with open(args.out_dir / f"{base}_tosconop_autodual_candidates.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    info["image"] = str(args.image)
    (args.out_dir / "extinction_selector_info.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    print(json.dumps(info), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
