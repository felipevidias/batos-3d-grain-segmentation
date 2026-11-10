"""Tile grid, frozen-method check and the per-tile worker (frozen v3 on a padded tile + centroid ownership).

Face rejection: frozen v3 rejects candidates touching the faces of the volume it is given. Run on the PADDED tile,
those faces are the outer halo faces (or the true volume faces at the volume border), so grains touching the
artificial core boundary are never rejected. No frozen file is edited.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import shutil
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import tifffile
import yaml

REPO = Path(__file__).resolve().parents[2]
MAIN_CHECKOUT = Path(r"C:\Users\USER\batos-3d-grain-segmentation")  # editable install target (subprocess stages)
METHOD_CONFIG = REPO / "configs" / "batos_paper_v3_scalefree.yaml"
ROUTE_FILES = {"BA-TOS_sel": ("e23c/watershed", "batos_selective"), "BA-TOS_ToS": ("e11/watershed", "tos_e11_min")}


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def check_freeze() -> dict:
    """Fail loudly unless every frozen file has its frozen SHA-256 in BOTH this worktree and the main checkout
    (frozen stages run as `python -I -m batos_v3...` subprocesses, which resolve to the main checkout)."""
    cfg = yaml.safe_load(METHOD_CONFIG.read_text(encoding="utf-8"))
    if cfg.get("status") != "frozen_for_test":
        raise RuntimeError(f"{METHOD_CONFIG} is not frozen_for_test")
    bad = []
    for k, v in cfg["code_sha256_at_freeze"].items():
        for root in (REPO, MAIN_CHECKOUT):
            if sha256(root / k) != v:
                bad.append(str(root / k))
    bad += [k for k, v in cfg["recovered_scripts_sha256"].items() if sha256(Path(v["path"])) != v["sha256"]]
    for rel in ("src/batos/batos_background.py", "src/batos/batos_markers.py"):  # worktree copies used by SCRIPTS
        key = {"src/batos/batos_background.py": "batos_background", "src/batos/batos_markers.py": "batos_markers"}[rel]
        if sha256(REPO / rel) != cfg["recovered_scripts_sha256"][key]["sha256"]:
            bad.append(str(REPO / rel))
    if bad:
        raise RuntimeError(f"code changed since the v3 freeze: {bad}")
    return cfg


@dataclass
class Tile:
    tile_id: str
    iz: int
    iy: int
    ix: int
    core: tuple  # ((z0, z1), (y0, y1), (x0, x1)) global, half-open
    pad: tuple   # padded box, clipped to the volume

    def core_slices(self):
        return tuple(slice(a, b) for a, b in self.core)

    def pad_slices(self):
        return tuple(slice(a, b) for a, b in self.pad)


def make_grid(shape, core: int, halo: int) -> list[Tile]:
    n = [int(np.ceil(s / core)) for s in shape]
    tiles = []
    for iz in range(n[0]):
        for iy in range(n[1]):
            for ix in range(n[2]):
                c = tuple((i * core, min((i + 1) * core, s)) for i, s in zip((iz, iy, ix), shape))
                p = tuple((max(0, a - halo), min(s, b + halo)) for (a, b), s in zip(c, shape))
                tiles.append(Tile(f"t{iz}_{iy}_{ix}", iz, iy, ix, c, p))
    return tiles


def select_tiles(tiles: list[Tile], shape, core: int, block: str | None, ids: list[str] | None) -> list[Tile]:
    if ids:
        sel = [t for t in tiles if t.tile_id in set(ids)]
    elif block == "all":
        sel = list(tiles)
    elif block and block.startswith("center"):
        k = int(block[len("center"):])  # e.g. center3 -> 3x3x3 tiles around the volume centre
        cen = [int(round(s / 2 / core - k / 2)) for s in shape]
        sel = [t for t in tiles if all(c <= i < c + k for c, i in zip(cen, (t.iz, t.iy, t.ix)))]
    else:
        raise ValueError("give --block center<k>|all or --tiles")
    if not sel:
        raise RuntimeError("no tile selected")
    return sel


def _log(msg: str) -> None:
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def instances_table(seg: np.ndarray, tile: Tile) -> pd.DataFrame:
    """Per instance of a padded-tile segmentation: voxels, global centroid / bbox, ownership by centroid."""
    off = np.array([a for a, _ in tile.pad])
    ids = np.flatnonzero(np.bincount(seg.ravel()))
    ids = ids[ids > 0]
    if len(ids) == 0:
        return pd.DataFrame(columns=["label", "voxels", "cz", "cy", "cx", "z0", "z1", "y0", "y1", "x0", "x1",
                                     "owned", "crosses_core"])
    vox = np.bincount(seg.ravel())[ids]
    com = np.asarray(ndi.center_of_mass(np.ones(seg.shape, np.uint8), seg, ids)) + off
    boxes = ndi.find_objects(seg)
    bb = np.array([[s.start for s in boxes[i - 1]] + [s.stop for s in boxes[i - 1]] for i in ids]) + np.r_[off, off]
    cfl = np.floor(com).astype(int)
    owned = np.all([(cfl[:, k] >= tile.core[k][0]) & (cfl[:, k] < tile.core[k][1]) for k in range(3)], axis=0)
    crosses = np.any([(bb[:, k] < tile.core[k][0]) | (bb[:, k + 3] > tile.core[k][1]) for k in range(3)], axis=0)
    return pd.DataFrame({"label": ids, "voxels": vox, "cz": com[:, 0], "cy": com[:, 1], "cx": com[:, 2],
                         "z0": bb[:, 0], "z1": bb[:, 3], "y0": bb[:, 1], "y1": bb[:, 4], "x0": bb[:, 2],
                         "x1": bb[:, 5], "owned": owned, "crosses_core": crosses})


def process_tile(tile: Tile, gray_path: Path, out_root: Path, halo: int, routes=("BA-TOS_sel",)) -> dict:
    """Run frozen v3 on the padded tile, keep owned instances, write kept labels + seed voxels, prune."""
    from batos_v3.pipeline import StageError, run_crop  # imported here: heavy, and only needed in workers

    out = out_root / "tiles" / tile.tile_id
    rec_path = out / "tile.json"
    if rec_path.exists():
        return json.loads(rec_path.read_text(encoding="utf-8"))
    if out.exists():
        shutil.rmtree(out)  # partial run of this driver (our own output)
    out.mkdir(parents=True)
    t0 = time.time()
    gray = tifffile.memmap(gray_path, mode="r")
    sub = np.ascontiguousarray(gray[tile.pad_slices()])
    del gray
    gtif = out / "gray.tif"
    tifffile.imwrite(gtif, sub)
    shape = sub.shape
    del sub
    rec = {"tile_id": tile.tile_id, **{k: v for k, v in asdict(tile).items() if k != "tile_id"},
           "padded_shape": list(shape), "halo": halo, "status": None}
    t1 = time.time()
    try:
        r = run_crop(tile.tile_id, gtif, out / "v3", "v3")
        rec["status"] = "EMPTY" if r.get("empty_result") else "DONE"
        g = r.get("gate") or {}
        rec |= {"r_granulo": r.get("r_granulo"), "s": r.get("s"), "size_filter_own": r.get("size_filter"),
                "openrec_radius": r.get("openrec_radius"), "openrec_impl": r.get("openrec_implementation"),
                "stage03_thr": r.get("stage03_dark_threshold"), "specimen_fraction": r.get("specimen_fraction"),
                "gate_triggered": g.get("gate_triggered"), "mv": g.get("region_voxels_median"),
                "qv": g.get("region_voxels_q75"), "tos_markers": r.get("tos_markers"),
                "v3_seconds": r.get("seconds_total"), "peak_memory_mb": r.get("peak_memory_mb"),
                "git": (r.get("git") or {}).get("commit"), "steps": r.get("steps")}
    except StageError as e:  # recorded as an empty tile (no tuning, no retry), as in run_test_v3
        rec |= {"status": "FAILED", "error": str(e), "v3_seconds": round(time.time() - t1, 1)}
        run_json = out / "v3" / "v3_run.json"
        steps_log = out / "v3" / "commands.log"
        rec["log_tail"] = steps_log.read_text(encoding="utf-8", errors="replace")[-1500:] if steps_log.exists() else ""
        if run_json.exists():
            rec["v3_run"] = json.loads(run_json.read_text(encoding="utf-8"))
    t2 = time.time()
    for route in routes:
        if rec["status"] != "DONE":
            rec[f"{route}_instances"] = 0
            rec[f"{route}_owned"] = 0
            continue
        ws, lab = ROUTE_FILES[route]
        seg = tifffile.imread(out / "v3" / ws / "segmentations" / f"{lab}_distance_c1.tif").astype(np.uint32)
        mk = tifffile.imread(out / "v3" / ws / "normalized_markers" / f"{lab}.tif")
        if seg.shape != shape or mk.shape != shape:
            raise RuntimeError(f"{tile.tile_id}: v3 output shape {seg.shape}/{mk.shape} != tile {shape}")
        tab = instances_table(seg, tile)
        owned = tab.loc[tab.owned, "label"].to_numpy()
        kept = np.where(np.isin(seg, owned), seg, 0).astype(np.uint32)
        tifffile.imwrite(out / f"kept_{route}.tif", kept, compression="zlib")
        tifffile.imwrite(out / f"full_{route}.tif", seg, compression="zlib")  # forensics: un-owned instances too
        mz = np.argwhere((mk > 0) & (kept > 0))
        np.savez_compressed(out / f"seeds_{route}.npz",
                            zyx=(mz + np.array([a for a, _ in tile.pad])).astype(np.int32),
                            label=kept[tuple(mz.T)].astype(np.uint32))
        tab.to_csv(out / f"instances_{route}.csv", index=False)
        rec[f"{route}_instances"] = int(len(tab))
        rec[f"{route}_owned"] = int(len(owned))
        rec[f"{route}_owned_crossing_core"] = int(tab.loc[tab.owned, "crosses_core"].sum())
    # prune: every v3 TIFF and the gray tile go; kept/full segmentations (zlib), seeds and JSON/CSV/logs stay
    for p in (out / "v3").rglob("*.tif"):
        p.unlink()
    gtif.unlink()
    rec |= {"post_seconds": round(time.time() - t2, 1), "seconds_total": round(time.time() - t0, 1),
            "finished": dt.datetime.now().astimezone().isoformat(timespec="seconds")}
    rec_path.write_text(json.dumps(rec, indent=1, default=str), encoding="utf-8")
    _log(f"{tile.tile_id} {rec['status']} {rec['seconds_total']:.0f}s r {rec.get('r_granulo')} gate "
         f"{rec.get('gate_triggered')} peak {rec.get('peak_memory_mb')} MB owned {rec.get('BA-TOS_sel_owned')}"
         f"/{rec.get('BA-TOS_sel_instances')}")
    return rec


if __name__ == "__main__":  # pragma: no cover
    sys.exit("use scripts/run_fullvolume.py")
