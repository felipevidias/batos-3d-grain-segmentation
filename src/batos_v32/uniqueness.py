"""Per-grain marker uniqueness rules applied to the e11 ToS-MIN marker list (Phase 5h). GT is never read.

Input: the e11 builder output (markers.csv: one row per point marker, column `priority`, `source_row_index`,
`projected_z/y/x`, `node`; markers_uint16.tif labelled 1..n in csv row order). Output: the same two files restricted
to the kept rows (same columns, same order, relabelled 1..m) so every downstream stage (EDT-Otsu count match,
watershed, E23c gate) consumes them unchanged.

Rules (one per hypothesis; s = r_granulo / R_REF):
  none : identity (must reproduce v3.1 bit-for-bit).
  H1   : one marker per basin of the EDT watershed of the Otsu mask. Mask = the e11 watershed/control Otsu mask
         (clean_bright_otsu_mask of the patched e11 control script, min component 500 s^3, closing round(s)), basin
         seeds = the EDT-Otsu control's own markers before count matching (EDT local maxima, greedy NMS 6 s), basins =
         skimage watershed(-EDT, seeds, mask). Per basin keep the highest-priority marker (ties: lower
         source_row_index), as the e11 builder does per internal component. Markers outside the mask (basin 0) are kept.
  H2   : greedy NMS on the marker list, radius 0.8 r_granulo, priority order (the e11 builder's greedy_nms).
  H3   : merge markers whose exclusive ToS supports (leaf labels of the e10 propagation, batos_v32.supports) touch
         (6-adjacency) and whose combined support fits one grain volume V_g = 4/3 pi r_granulo^3. Agglomerative:
         touching pairs in increasing order of combined size, union-find, accepted while the union is <= V_g. Each
         merged group keeps its highest-priority marker.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

PATCHED_CONTROL = Path(__file__).resolve().parents[1] / "batos_v3" / "patched" / "e11_build_control_markers_v3.py"
H2_FACTOR = 0.8


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def priority_order(m: pd.DataFrame) -> np.ndarray:
    return np.lexsort((m["source_row_index"].to_numpy(np.int64), -m["priority"].to_numpy(float)))


def points(m: pd.DataFrame) -> np.ndarray:
    return m[["projected_z", "projected_y", "projected_x"]].to_numpy(np.int64)


def h1_keep(m: pd.DataFrame, flat: np.ndarray, min_component: int, closing: int, nms: float) -> tuple[np.ndarray, dict]:
    from skimage.segmentation import watershed
    ctl = _load(PATCHED_CONTROL, "e11_control_v3")
    mask, thr = ctl.clean_bright_otsu_mask(ctl.normalize_image(flat), min_component, closing)
    seeds = ctl.greedy_nms(ctl.otsu_edt_candidates(mask), radius=float(nms))
    edt = ndi.distance_transform_edt(mask)
    sv = np.zeros(mask.shape, np.int32)
    p = seeds[["z", "y", "x"]].to_numpy(np.int64)
    sv[tuple(p.T)] = np.arange(1, len(p) + 1, dtype=np.int32)
    basins = watershed(-edt, sv, mask=mask)
    del edt, sv
    pts = points(m)
    b = basins[tuple(pts.T)].astype(np.int64)
    keep = np.zeros(len(m), bool)
    seen: set[int] = set()
    for i in priority_order(m):
        if b[i] == 0:
            keep[i] = True
        elif b[i] not in seen:
            seen.add(int(b[i]))
            keep[i] = True
    occ = np.bincount(b[b > 0], minlength=len(p) + 1)[1:]
    return keep, {"otsu_threshold": float(thr), "basins": int(len(p)), "markers_outside_mask": int((b == 0).sum()),
                  "basins_with_markers": int((occ > 0).sum()), "basins_with_multiple_markers": int((occ > 1).sum())}


def h2_keep(m: pd.DataFrame, radius: float) -> tuple[np.ndarray, dict]:
    pts = points(m).astype(float)
    r2 = radius * radius
    keep = np.zeros(len(m), bool)
    kept: list[int] = []
    for i in priority_order(m):  # O(n * kept); n <= a few thousand
        if kept:
            d = pts[kept] - pts[i]
            if (np.einsum("ij,ij->i", d, d) < r2).any():
                continue
        kept.append(int(i))
        keep[i] = True
    return keep, {"nms_radius": radius}


def touching_pairs(leaf: np.ndarray, labels: np.ndarray) -> set[tuple[int, int]]:
    """6-adjacent pairs (a < b) of leaf labels, restricted to `labels`."""
    lut = np.zeros(int(leaf.max()) + 1, bool)
    lut[labels] = True
    pairs: set[tuple[int, int]] = set()
    for ax in range(3):
        a = np.moveaxis(leaf, ax, 0)
        u, v = a[:-1], a[1:]
        sel = (u != v) & lut[u] & lut[v]
        if sel.any():
            uv = np.stack([u[sel], v[sel]], 1)
            uv.sort(axis=1)
            uv = np.unique(uv, axis=0)
            pairs.update(map(tuple, uv.tolist()))
    return pairs


def h3_keep(m: pd.DataFrame, diag: pd.DataFrame, leaf: np.ndarray, v_grain: float) -> tuple[np.ndarray, dict]:
    lab = diag.set_index("node")["candidate_label"]
    vox = diag.set_index("node")["region_voxels"]
    nodes = m["node"].to_numpy(np.int64)
    cl = lab.loc[nodes].to_numpy(np.int64)
    size = vox.loc[nodes].to_numpy(np.int64)
    idx = {int(c): i for i, c in enumerate(cl)}
    pairs = touching_pairs(leaf, cl)
    edges = sorted(((size[idx[a]] + size[idx[b]], idx[a], idx[b]) for a, b in pairs))
    parent = list(range(len(m)))
    tot = size.astype(np.int64).copy()

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    merges = 0
    for _, i, j in edges:
        ri, rj = find(i), find(j)
        if ri == rj or tot[ri] + tot[rj] > v_grain:
            continue
        parent[rj] = ri
        tot[ri] += tot[rj]
        merges += 1
    keep = np.zeros(len(m), bool)
    seen: set[int] = set()
    for i in priority_order(m):
        r = find(int(i))
        if r not in seen:
            seen.add(r)
            keep[i] = True
    return keep, {"v_grain": v_grain, "touching_pairs": len(pairs), "merges": merges,
                  "support_voxels_median": float(np.median(size)) if len(size) else None}


def write(m: pd.DataFrame, keep: np.ndarray, shape, out: Path, info: dict) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    k = m.loc[keep].reset_index(drop=True)
    k.to_csv(out / "markers.csv", index=False)
    vol = np.zeros(shape, np.uint16 if len(k) <= 65535 else np.uint32)
    p = points(k)
    if len(k):
        vol[tuple(p.T)] = np.arange(1, len(k) + 1)
    tifffile.imwrite(out / "markers_uint16.tif", vol, compression="zlib")
    info = {"markers_in": int(len(m)), "markers_out": int(len(k)), "gt_used": False} | info
    (out / "uniqueness.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    return info


def apply(hyp: str, tos_dir: Path, out: Path, *, r_granulo: float, s: float, flat: Path | None = None,
          diag: Path | None = None, leaf: Path | None = None, min_component: int = 500, closing: int = 1,
          edt_nms: float = 6.0) -> dict:
    m = pd.read_csv(tos_dir / "markers.csv")
    shape = tifffile.imread(tos_dir / "markers_uint16.tif").shape
    info: dict = {"hypothesis": hyp, "r_granulo": r_granulo, "s": s}
    if hyp == "none":
        keep = np.ones(len(m), bool)
    elif hyp == "H1":
        keep, extra = h1_keep(m, tifffile.imread(flat), min_component, closing, edt_nms)
        info |= extra | {"min_component": min_component, "closing": closing, "edt_nms": edt_nms}
    elif hyp == "H2":
        keep, extra = h2_keep(m, H2_FACTOR * r_granulo)
        info |= extra
    elif hyp == "H3":
        keep, extra = h3_keep(m, pd.read_csv(diag, usecols=["node", "candidate_label", "region_voxels"]),
                              tifffile.imread(leaf), 4.0 / 3.0 * math.pi * r_granulo ** 3)
        info |= extra
    else:
        raise ValueError(hyp)
    return write(m, keep, shape, out, info)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--hypothesis", choices=["none", "H1", "H2", "H3"], required=True)
    ap.add_argument("--tos-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--r-granulo", type=float, required=True)
    ap.add_argument("--s", type=float, required=True)
    ap.add_argument("--flat", type=Path)
    ap.add_argument("--diag", type=Path)
    ap.add_argument("--leaf", type=Path)
    ap.add_argument("--min-component", type=int, default=500)
    ap.add_argument("--closing", type=int, default=1)
    ap.add_argument("--edt-nms", type=float, default=6.0)
    a = ap.parse_args(argv)
    info = apply(a.hypothesis, a.tos_dir, a.out, r_granulo=a.r_granulo, s=a.s, flat=a.flat, diag=a.diag,
                 leaf=a.leaf, min_component=a.min_component, closing=a.closing, edt_nms=a.edt_nms)
    print(json.dumps(info))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
