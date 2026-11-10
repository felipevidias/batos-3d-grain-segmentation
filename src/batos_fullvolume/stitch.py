"""Deterministic sequential stitch of the per-tile kept instances into one uint32 label volume (memmapped .npy).

Tiles are merged in index order (z, y, x). For tile T, every kept instance A is compared with the global labels
already written (instances of earlier tiles) under A's voxels:
  * one-to-one matching between A's and existing instances, greedy by descending IoU, pairs with IoU >= match_iou
    (0.5: a pair above it is unique on both sides, so the matching is one-to-one by construction);
    matched pair -> the larger instance is kept (ties -> the existing one); the smaller is deleted entirely;
  * unmatched overlap ("contested voxels"): a voxel inside T's core goes to A (T is authoritative in its core);
    otherwise it stays with the existing instance.
Global ids are assigned in stitch order (1, 2, ...). Only the domain box (union of the selected padded tiles) is
materialised; one bbox-sized block per instance is in RAM at a time.
Outputs (out_dir): labels.npy (+ labels_meta.json with the global offset), instances.csv (one row per global
instance, provenance + reconciliation flags), dropped.csv, seeds.npz (global zyx + gid of every seed voxel of a
surviving instance; at evaluation the support of p is the subset still labelled p).

Optional orphan rescue (orphan_rescue=True; added 2026-11-11 after the smoke test showed seam grains owned by no
tile: each tile's prediction of a plane-straddling grain is biased toward its own halo, so both centroids fall in
the other tile's core). After the owned pass, every NON-owned tile instance (tile order, then local label) whose
centroid lies in the core of a selected tile and whose voxels are still mostly unlabelled (free fraction >=
rescue_min_free) is written into the FREE voxels only (never steals) with a new gid, flag rescued=1. The second
tile's copy of the same grain then finds its voxels taken and is skipped. Seeds: marker voxels of the instance
from seeds_all_<route>.npz when the tile has it; otherwise (tiles run before that file existed) the free instance
voxel nearest the centroid is used as a single proxy seed (seed_proxy=1, an approximation for seed-ignore only).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import tifffile

from batos_fullvolume.tiling import Tile


def domain_box(tiles: list[Tile]) -> tuple:
    return tuple((min(t.pad[k][0] for t in tiles), max(t.pad[k][1] for t in tiles)) for k in range(3))


def stitch(tiles: list[Tile], tile_root: Path, out_dir: Path, route: str = "BA-TOS_sel",
           match_iou: float = 0.5, orphan_rescue: bool = False, rescue_min_free: float = 0.5) -> dict:
    t_start = time.time()
    out_dir.mkdir(parents=True, exist_ok=True)
    box = domain_box(tiles)
    off = np.array([a for a, _ in box])
    shape = tuple(b - a for a, b in box)
    lab_path = out_dir / "labels.npy"
    L = np.lib.format.open_memmap(lab_path, mode="w+", dtype=np.uint32, shape=shape)
    reg: dict[int, dict] = {}       # gid -> live record
    dropped: list[dict] = []
    seeds_zyx, seeds_gid = [], []
    gid_next = 1
    for tile in sorted(tiles, key=lambda t: (t.iz, t.iy, t.ix)):
        trec = json.loads((tile_root / tile.tile_id / "tile.json").read_text(encoding="utf-8"))
        if trec["status"] != "DONE":
            continue
        kept = tifffile.imread(tile_root / tile.tile_id / f"kept_{route}.tif")
        inst = pd.read_csv(tile_root / tile.tile_id / f"instances_{route}.csv")
        inst = inst[inst.owned].set_index("label")
        sd = np.load(tile_root / tile.tile_id / f"seeds_{route}.npz")
        pad0 = np.array([a for a, _ in tile.pad])
        boxes = ndi.find_objects(kept)
        core_lo = np.array([a for a, _ in tile.core]) - off
        core_hi = np.array([b for _, b in tile.core]) - off
        # 1) overlaps of every kept instance with what is already written (before writing anything of T)
        blocks, pairs = {}, []
        for a in inst.index:
            sl_t = boxes[a - 1]
            g0 = pad0 + [s.start for s in sl_t] - off
            g1 = pad0 + [s.stop for s in sl_t] - off
            sl_g = tuple(slice(int(x), int(y)) for x, y in zip(g0, g1))
            mask = kept[sl_t] == a
            cur = np.asarray(L[sl_g])
            under = cur[mask]
            ids, cnt = np.unique(under[under > 0], return_counts=True)
            blocks[a] = (sl_g, mask, g0)
            na = int(mask.sum())
            for b, n in zip(ids.tolist(), cnt.tolist()):
                nb = reg[b]["voxels"]
                pairs.append((n / (na + nb - n), a, b, n))
        # 2) one-to-one greedy matching at IoU >= match_iou
        pairs.sort(key=lambda x: (-x[0], x[1], x[2]))
        used_a, used_b, match = set(), set(), {}
        for iou, a, b, n in pairs:
            if iou >= match_iou and a not in used_a and b not in used_b:
                used_a.add(a), used_b.add(b)
                match[a] = (b, iou)
        overlapped = {a for _, a, _, _ in pairs}
        # 3) write
        for a in inst.index:
            sl_g, mask, g0 = blocks[a]
            na = int(mask.sum())
            row = {"tile_id": tile.tile_id, "local_label": int(a), "tile_voxels": na,
                   "cz": float(inst.at[a, "cz"]), "cy": float(inst.at[a, "cy"]), "cx": float(inst.at[a, "cx"]),
                   "crosses_core": bool(inst.at[a, "crosses_core"]), "s_owner": trec.get("s"),
                   "own_filter": trec.get("size_filter_own"),
                   "matched_replaced": 0, "matched_iou": np.nan, "contested_won": 0, "contested_lost": 0,
                   "overlap_any": a in overlapped, "rescued": 0, "seed_proxy": 0}
            if a in match:
                b, iou = match[a]
                if na <= reg[b]["voxels"]:  # existing is at least as large -> keep it, drop A
                    dropped.append({**row, "reason": "matched_smaller", "kept_gid": b, "iou": iou})
                    reg[b]["matched_replaced"] += 1
                    reg[b]["matched_iou"] = iou
                    continue
                # A larger: delete B entirely
                bsl = reg[b]["slices"]
                blk = np.asarray(L[bsl])
                blk[blk == b] = 0
                L[bsl] = blk
                dropped.append({**{k: v for k, v in reg[b].items() if k != "slices"}, "reason": "matched_larger_new",
                                "kept_gid": gid_next, "iou": iou})
                del reg[b]
                row["matched_replaced"] = 1
                row["matched_iou"] = iou
            cur = np.asarray(L[sl_g])
            zz, yy, xx = np.nonzero(mask)
            coords = np.stack([zz, yy, xx], 1) + g0
            in_core = np.all((coords >= core_lo) & (coords < core_hi), axis=1)
            curv = cur[mask]
            take = (curv == 0) | ((curv > 0) & in_core)
            lost = curv[take & (curv > 0)]
            if len(lost):
                lb, ln = np.unique(lost, return_counts=True)
                for b, n in zip(lb.tolist(), ln.tolist()):
                    reg[b]["voxels"] -= n
                    reg[b]["contested_lost"] += n
            row["contested_won"] = int(len(lost))
            gid = gid_next
            gid_next += 1
            newv = curv.copy()
            newv[take] = gid
            cur[mask] = newv
            L[sl_g] = cur
            row["voxels"] = int(take.sum())
            reg[gid] = {"gid": gid, **row, "slices": sl_g}
            s = sd["label"] == a
            seeds_zyx.append(sd["zyx"][s])
            seeds_gid.append(np.full(int(s.sum()), gid, np.uint32))
    L.flush()
    n_rescue_cand = 0
    if orphan_rescue:
        sel_idx = {(t.iz, t.iy, t.ix) for t in tiles}
        core = tiles[0].core[0][1] - tiles[0].core[0][0]
        for tile in sorted(tiles, key=lambda t: (t.iz, t.iy, t.ix)):
            d = tile_root / tile.tile_id
            trec = json.loads((d / "tile.json").read_text(encoding="utf-8"))
            if trec["status"] != "DONE":
                continue
            full = tifffile.imread(d / f"full_{route}.tif")
            tab = pd.read_csv(d / f"instances_{route}.csv").set_index("label")
            sda = np.load(d / f"seeds_all_{route}.npz") if (d / f"seeds_all_{route}.npz").exists() else None
            pad0 = np.array([a for a, _ in tile.pad])
            boxes = ndi.find_objects(full)
            for a in tab.index[~tab.owned.to_numpy(bool)]:
                c = tab.loc[a, ["cz", "cy", "cx"]].to_numpy(float)
                if tuple(int(v) // core for v in np.floor(c)) not in sel_idx:
                    continue  # centroid outside the selected cores: a grain of an unselected tile
                sl_t = boxes[a - 1]
                g0 = pad0 + [s_.start for s_ in sl_t] - off
                g1 = pad0 + [s_.stop for s_ in sl_t] - off
                if np.any(g0 < 0) or np.any(g1 > np.array(shape)):
                    continue
                n_rescue_cand += 1
                sl_g = tuple(slice(int(x), int(y)) for x, y in zip(g0, g1))
                mask = full[sl_t] == a
                cur = np.asarray(L[sl_g])
                free = mask & (cur == 0)
                na, nf = int(mask.sum()), int(free.sum())
                if nf < rescue_min_free * na:
                    continue
                gid = gid_next
                gid_next += 1
                cur[free] = gid
                L[sl_g] = cur
                reg[gid] = {"gid": gid, "tile_id": tile.tile_id, "local_label": int(a), "tile_voxels": na,
                            "cz": float(c[0]), "cy": float(c[1]), "cx": float(c[2]),
                            "crosses_core": bool(tab.loc[a, "crosses_core"]), "s_owner": trec.get("s"),
                            "own_filter": trec.get("size_filter_own"), "matched_replaced": 0, "matched_iou": np.nan,
                            "contested_won": 0, "contested_lost": 0, "overlap_any": nf < na, "rescued": 1,
                            "seed_proxy": 0, "voxels": nf, "slices": sl_g}
                z_ = np.zeros((0, 3), np.int32)
                if sda is not None:
                    z_ = sda["zyx"][sda["label"] == a]
                    if len(z_):
                        z_ = z_[np.asarray(L[tuple((z_ - off).T)]) == gid]
                if len(z_) == 0:  # proxy seed: the free instance voxel nearest the centroid
                    fz = np.argwhere(free) + g0 + off
                    z_ = fz[[int(np.argmin(((fz - c) ** 2).sum(1)))]]
                    reg[gid]["seed_proxy"] = 1
                seeds_zyx.append(z_.astype(np.int32))
                seeds_gid.append(np.full(len(z_), gid, np.uint32))
        L.flush()
    # orphan diagnostic: every NON-owned tile instance vs the final labels (best IoU). An instance whose grain
    # no tile owns (each tile's centroid falls in the other's core) has best IoU ~0: a seam loss.
    orphans = []
    for tile in sorted(tiles, key=lambda t: (t.iz, t.iy, t.ix)):
        d = tile_root / tile.tile_id
        if json.loads((d / "tile.json").read_text(encoding="utf-8"))["status"] != "DONE":
            continue
        full = tifffile.imread(d / f"full_{route}.tif")
        tab = pd.read_csv(d / f"instances_{route}.csv")
        pad0 = np.array([a for a, _ in tile.pad])
        boxes = ndi.find_objects(full)
        for a in tab.loc[~tab.owned, "label"]:
            sl_t = boxes[a - 1]
            g0 = pad0 + [s.start for s in sl_t] - off
            g1 = pad0 + [s.stop for s in sl_t] - off
            if np.any(g0 < 0) or np.any(g1 > np.array(shape)):
                continue  # instance leaves the stitched box (block border): not a seam case
            mask = full[sl_t] == a
            cur = np.asarray(L[tuple(slice(int(x), int(y)) for x, y in zip(g0, g1))])[mask]
            ids, cnt = np.unique(cur[cur > 0], return_counts=True)
            best = 0.0
            for b, n in zip(ids.tolist(), cnt.tolist()):
                if b in reg:
                    best = max(best, n / (mask.sum() + reg[b]["voxels"] - n))
            orphans.append({"tile_id": tile.tile_id, "local_label": int(a), "voxels": int(mask.sum()),
                            "best_iou_final": best})
    orph = pd.DataFrame(orphans)
    orph.to_csv(out_dir / "nonowned_vs_final.csv", index=False)
    del L
    inst_df = pd.DataFrame([{k: v for k, v in r.items() if k != "slices"} for r in reg.values()])
    inst_df.to_csv(out_dir / "instances.csv", index=False)
    pd.DataFrame(dropped).to_csv(out_dir / "dropped.csv", index=False)
    zyx = np.concatenate(seeds_zyx) if seeds_zyx else np.zeros((0, 3), np.int32)
    gids = np.concatenate(seeds_gid) if seeds_gid else np.zeros(0, np.uint32)
    alive = np.isin(gids, inst_df["gid"].to_numpy()) if len(inst_df) else np.zeros(0, bool)
    np.savez_compressed(out_dir / "seeds.npz", zyx=zyx[alive], gid=gids[alive])
    meta = {"route": route, "match_iou": match_iou, "contested_voxels": "core_authority",
            "box": [list(b) for b in box], "offset": off.tolist(), "shape": list(shape), "dtype": "uint32",
            "labels": str(lab_path), "n_instances": int(len(inst_df)), "n_dropped": len(dropped),
            "max_gid": gid_next - 1,
            "orphan_rescue": orphan_rescue, "rescue_min_free": rescue_min_free if orphan_rescue else None,
            "rescue_candidates": n_rescue_cand,
            "n_rescued": int(sum(r["rescued"] for r in reg.values())),
            "n_rescued_seed_proxy": int(sum(r["seed_proxy"] for r in reg.values())),
            "nonowned_checked": int(len(orph)),
            "nonowned_best_iou_lt050": int((orph.best_iou_final < 0.5).sum()) if len(orph) else 0,
            "nonowned_best_iou_lt050_ge1000vox": int(((orph.best_iou_final < 0.5) & (orph.voxels >= 1000)).sum())
            if len(orph) else 0,
            "tiles": [t.tile_id for t in tiles],
            "seconds": round(time.time() - t_start, 1)}
    (out_dir / "labels_meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    return meta
