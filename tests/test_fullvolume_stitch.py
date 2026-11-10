"""Synthetic checks of batos_fullvolume.stitch (matching, keep-larger, core authority, seeds) and of the
ownership rule in batos_fullvolume.tiling.instances_table."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from batos_fullvolume.stitch import stitch  # noqa: E402
from batos_fullvolume.tiling import instances_table, make_grid  # noqa: E402


def _write_tile(root: Path, tile, seg: np.ndarray):
    d = root / tile.tile_id
    d.mkdir(parents=True)
    tab = instances_table(seg, tile)
    owned = tab.loc[tab.owned, "label"].to_numpy()
    kept = np.where(np.isin(seg, owned), seg, 0).astype(np.uint32)
    tifffile.imwrite(d / "kept_BA-TOS_sel.tif", kept)
    tifffile.imwrite(d / "full_BA-TOS_sel.tif", seg)
    tab.to_csv(d / "instances_BA-TOS_sel.csv", index=False)
    z = np.argwhere(kept > 0)  # every voxel a seed (simple)
    np.savez(d / "seeds_BA-TOS_sel.npz", zyx=(z + [a for a, _ in tile.pad]).astype(np.int32),
             label=kept[tuple(z.T)].astype(np.uint32))
    (d / "tile.json").write_text(json.dumps({"status": "DONE", "s": 1.0, "size_filter_own": 500}))
    return tab


def test_stitch(tmp_path):
    shape, core, halo = (20, 10, 10), 10, 4
    tiles = make_grid(shape, core, halo)  # t0_0_0 core z 0-10, t1_0_0 core z 10-20
    t0, t1 = tiles
    assert t0.pad == ((0, 14), (0, 10), (0, 10)) and t1.pad == ((6, 20), (0, 10), (0, 10))
    # global truth: grain A z 7..13 (centroid 9.5 -> tile 0), grain B z 15..19 (tile 1)
    s0 = np.zeros((14, 10, 10), np.uint32)
    s0[7:14, 2:6, 2:6] = 5          # tile 0 sees A as z 7..13 (centroid 10 -> owned by tile 1!)
    s0[1:4, 1:4, 1:4] = 6           # interior grain of tile 0
    s1 = np.zeros((14, 10, 10), np.uint32)
    s1[1:7, 2:6, 2:6] = 3           # tile 1 sees A as z 7..12 (centroid 9.5 -> tile 0)
    s1[9:13, 6:9, 6:9] = 4          # grain B, z 15..18, owned by tile 1
    tab0 = _write_tile(tmp_path / "tiles", t0, s0)
    tab1 = _write_tile(tmp_path / "tiles", t1, s1)
    assert tab0.set_index("label").loc[5, "owned"] == False  # noqa: E712  centroid z 10.0 -> core of tile 1
    assert tab1.set_index("label").loc[3, "owned"] == False  # noqa: E712  centroid z 9.5 -> core of tile 0
    # -> A is lost by both tiles (ownership disagreement); check it is reported as such: no instance for A
    meta = stitch(tiles, tmp_path / "tiles", tmp_path / "st")
    L = np.load(tmp_path / "st" / "labels.npy")
    assert meta["n_instances"] == 2 and set(np.unique(L)) == {0, 1, 2}
    assert meta["nonowned_best_iou_lt050"] == 2  # A seen by both tiles, owned by none -> two orphan records


def test_stitch_match_and_core_authority(tmp_path):
    shape, core, halo = (20, 10, 10), 10, 4
    t0, t1 = make_grid(shape, core, halo)
    s0 = np.zeros((14, 10, 10), np.uint32)
    s0[6:13, 2:6, 2:6] = 5          # global z 6..12, centroid 9 -> owned by tile 0 (112 voxels)
    s1 = np.zeros((14, 10, 10), np.uint32)
    s1[0:5, 2:6, 2:6] = 3           # global z 6..10, centroid 8 -> tile 0 -> not owned
    s1[4:8, 2:6, 2:7] = 7           # global z 10..13, x 2..6: centroid 11.5 -> owned by tile 1 (80 voxels)
    _write_tile(tmp_path / "tiles", t0, s0)
    _write_tile(tmp_path / "tiles", t1, s1)
    stitch([t0, t1], tmp_path / "tiles", tmp_path / "st")
    L = np.load(tmp_path / "st" / "labels.npy")
    inst = pd.read_csv(tmp_path / "st" / "instances.csv").set_index("gid")
    # overlap z 10..12, x 2..5: 3*4*4 = 48 voxels; IoU = 48 / (112 + 80 - 48) = 0.33 < 0.5 -> not matched;
    # contested voxels z 10..12 lie in tile 1's core -> go to instance 2 (tile 1)
    assert (L[10:13, 2:6, 2:6] == 2).all() and (L[6:10, 2:6, 2:6] == 1).all()
    assert inst.loc[1, "contested_lost"] == 48 and inst.loc[2, "contested_won"] == 48
    assert inst.loc[1, "voxels"] == 112 - 48
    seeds = np.load(tmp_path / "st" / "seeds.npz")
    assert set(np.unique(seeds["gid"])) == {1, 2}


def test_stitch_matched_keep_larger(tmp_path):
    shape, core, halo = (20, 10, 10), 10, 4
    t0, t1 = make_grid(shape, core, halo)
    s0 = np.zeros((14, 10, 10), np.uint32)
    s0[5:13, 2:6, 2:6] = 5          # global z 5..12 (128 voxels), centroid 8.5 -> tile 0
    s1 = np.zeros((14, 10, 10), np.uint32)
    s1[4:9, 2:6, 2:6] = 7           # global z 10..14, 80 vox, centroid 12 -> tile 1
    s1[4:9, 6:7, 2:6] = 7           # + 20 voxels -> 100 voxels; overlap with A = z 10..12 x 4x4 = 48
    s1[0:4, 2:6, 2:6] = 7           # extend down to global z 6..9 -> 64 more: 164 voxels, overlap 112
    _write_tile(tmp_path / "tiles", t0, s0)
    tab1 = _write_tile(tmp_path / "tiles", t1, s1)
    assert tab1.set_index("label").loc[7, "owned"]  # centroid z ~10.1
    stitch([t0, t1], tmp_path / "tiles", tmp_path / "st")
    L = np.load(tmp_path / "st" / "labels.npy")
    inst = pd.read_csv(tmp_path / "st" / "instances.csv")
    drop = pd.read_csv(tmp_path / "st" / "dropped.csv")
    # IoU = 112 / (128 + 164 - 112) = 0.62 >= 0.5 -> matched; tile-1 instance larger -> replaces gid 1
    assert len(inst) == 1 and inst.gid.iloc[0] == 2 and inst.matched_replaced.iloc[0] == 1
    assert drop.reason.tolist() == ["matched_larger_new"]
    assert not (L == 1).any() and (L == 2).sum() == 164
