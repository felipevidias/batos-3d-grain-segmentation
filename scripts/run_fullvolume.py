#!/usr/bin/env python3
"""Phase 5g: full-volume tiled BA-TOS v3 (frozen method; configs/fullvolume_v3.yaml).

Subcommands (run in this order; each is restartable and skips finished work):
  scale     r_granulo once on the central 200^3 crop (GT-free, batos_v3.scale.estimate) -> halo = ceil(2 r)
  run       frozen v3 on every selected padded tile (<= --workers threads, each runs the frozen stage subprocesses;
            a new tile starts only while the estimated concurrent peak stays under --mem-budget-gb);
            writes tiles/<id>/{tile.json, kept_*.tif, full_*.tif, seeds_*.npz, instances_*.csv}, tiles_manifest.csv,
            resources.csv (process-tree RSS sampled every 2 s)
  stitch    deterministic merge -> stitch_<route>/labels.npy (+ instances.csv, dropped.csv, seeds.npz)
  evaluate  streaming evaluation against the full GT (valid = ROI labelsInside, >= 1000 voxels, inside the evaluated
            core box) -> eval_<route>/{metrics.csv, metrics.json, per_gt.csv, per_pred.csv, per_tile.csv, seam.json}
  report    heat map + REPORT.md
Example: python -I scripts/run_fullvolume.py run --name smoke_center3 --block center3 --workers 3
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import platform
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
import tifffile
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_fullvolume.tiling import (Tile, check_freeze, make_grid, process_tile,  # noqa: E402
                                     select_tiles, specimen_selection)

CFG = yaml.safe_load((REPO / "configs" / "fullvolume_v3.yaml").read_text(encoding="utf-8"))
GRAY = Path(CFG["volume"]["gray"])
GT = Path(CFG["volume"]["gt"])
ROI = Path(CFG["volume"]["roi_labels_inside"])
ROI_OUT = ROI.with_name(ROI.name.replace("Inside", "Outside"))
ROOT = Path(CFG["out_root"])
CORE = int(CFG["tiling"]["core"])


def log(msg: str) -> None:
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def git_commit() -> dict:
    def g(*a):
        return subprocess.run(["git", *a], cwd=REPO, capture_output=True, text=True).stdout.strip()
    return {"commit": g("rev-parse", "HEAD"), "branch": g("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty": bool(g("status", "--porcelain", "--untracked-files=no"))}


def vol_shape() -> tuple:
    with tifffile.TiffFile(GRAY) as t:
        return tuple(t.series[0].shape)


# ------------------------------------------------------------------------------------------------ scale
def cmd_scale(a) -> int:
    from batos_v3.scale import estimate
    out = ROOT / "scale_central.json"
    if out.exists() and not a.force:
        log(f"exists: {out}")
        print(out.read_text(encoding="utf-8"))
        return 0
    shape = vol_shape()
    c0 = [s // 2 - CORE // 2 for s in shape]
    box = [[c, c + CORE] for c in c0]
    gray = tifffile.memmap(GRAY, mode="r")
    crop = np.ascontiguousarray(gray[tuple(slice(*b) for b in box)])
    t0 = time.time()
    res, spec = estimate(crop)
    r = float(res["r_granulo"])
    halo = int(math.ceil(CFG["tiling"]["halo_factor"] * r))
    rec = {"volume": str(GRAY), "shape": list(shape), "central_box_zyx": box, "r_granulo": r, "r0": res["r0"],
           "specimen_fraction": res["specimen_fraction"], "halo_factor": CFG["tiling"]["halo_factor"], "halo": halo,
           "seconds": round(time.time() - t0, 1), "gt_used": False, "git": git_commit(),
           "timestamp": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "pattern_spectrum": res["pattern_spectrum"]}
    ROOT.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1), encoding="utf-8")
    log(f"r_granulo {r} (r0 {res['r0']}, specimen {res['specimen_fraction']:.3f}) -> halo {halo}")
    return 0


def load_scale() -> dict:
    p = ROOT / "scale_central.json"
    if not p.exists():
        raise FileNotFoundError(f"{p} missing: run `scale` first")
    return json.loads(p.read_text(encoding="utf-8"))


def run_dir(name: str) -> Path:
    return ROOT / name


def load_tiles(name: str) -> tuple[list[Tile], dict]:
    global CORE
    cfg = json.loads((run_dir(name) / "run_config.json").read_text(encoding="utf-8"))
    CORE = int(cfg["core"])
    tiles = [Tile(t["tile_id"], t["iz"], t["iy"], t["ix"], tuple(map(tuple, t["core"])), tuple(map(tuple, t["pad"])))
             for t in cfg["tiles"]]
    return tiles, cfg


# ------------------------------------------------------------------------------------------------ run
class Sampler(threading.Thread):
    def __init__(self, path: Path, running: callable):
        super().__init__(daemon=True)
        self.path, self.running, self.stop_ev, self.peak_mb = path, running, threading.Event(), 0.0

    def run(self):
        me = psutil.Process()
        new = not self.path.exists()
        with open(self.path, "a", encoding="utf-8") as f:
            if new:
                f.write("time,tree_rss_mb,n_running,system_available_mb\n")
            while not self.stop_ev.is_set():
                tot = 0
                for p in [me, *me.children(recursive=True)]:
                    try:
                        tot += p.memory_info().rss
                    except psutil.Error:
                        pass
                mb = tot / 2 ** 20
                self.peak_mb = max(self.peak_mb, mb)
                f.write(f"{dt.datetime.now():%Y-%m-%dT%H:%M:%S},{mb:.0f},{self.running()},"
                        f"{psutil.virtual_memory().available / 2 ** 20:.0f}\n")
                f.flush()
                self.stop_ev.wait(2.0)


class StageGate:
    """Memory admission per frozen-stage SUBPROCESS (scheduling only: command lines and frozen files unchanged).
    Every stage reserves est(stage) = 1.1 x the largest peak measured so far for that stage script (default
    `default_mb` until measured) and waits until the reservations of all running stages fit in `budget_mb`.
    Installed by replacing batos_v3.pipeline.run_tracked (the frozen launcher) with a gated wrapper in this
    process."""

    def __init__(self, budget_mb: float, default_mb: float, seen: dict | None = None):
        self.budget, self.default = budget_mb, default_mb
        self.used, self.cv = 0.0, threading.Condition()
        self.peak = dict(seen or {})

    @staticmethod
    def key(cmd: list) -> str:
        cmd = [str(c) for c in cmd]
        if "-m" in cmd:
            return cmd[cmd.index("-m") + 1]
        py = [c for c in cmd if c.endswith(".py")]
        return Path(py[0]).name if py else cmd[0]

    def est(self, k: str) -> float:
        return 1.1 * self.peak[k] if k in self.peak else self.default

    def install(self):
        import batos_v3.pipeline as pl
        orig = pl.run_tracked
        if getattr(orig, "_gated", False):
            return

        def gated(cmd, log_path):
            k = self.key(cmd)
            with self.cv:
                need = min(self.est(k), self.budget)
                while self.used + need > self.budget and self.used > 0:
                    self.cv.wait()
                self.used += need
            try:
                rc, peak = orig(cmd, log_path)
            finally:
                with self.cv:
                    self.used -= need
                    self.cv.notify_all()
            with self.cv:
                self.peak[k] = max(self.peak.get(k, 0.0), peak)
            return rc, peak
        gated._gated = True
        pl.run_tracked = gated


def install_scale_rule(central: dict, log_note) -> None:
    """Full-volume rule (decided 2026-11-10): scale is a specimen property, estimated ONCE per scan.

    The frozen scale stage still runs on every padded tile (it also writes the tile's specimen.tif, which
    v3 uses for its stage-03 threshold). Right after it, scale.json is rewritten with the per-scan
    r_granulo of the central 200^3 crop and openrec_radius = max(8, round(1.5 r)) (the batos_v3.scale
    formula); the tile's own estimate is kept as scale_tile_estimate.json. No frozen file is edited:
    batos_v3.pipeline.run_tracked is wrapped in this process, like the stage gate."""
    import batos_v3.pipeline as pl
    orig = pl.run_tracked
    r = float(central["r_granulo"])
    r_out = int(r) if r.is_integer() else r

    def wrapped(cmd, log_path):
        rc, peak = orig(cmd, log_path)
        c = [str(x) for x in cmd]
        if rc == 0 and "-m" in c and c[c.index("-m") + 1] == "batos_v3.scale":
            p = Path(c[c.index("--out-json") + 1])
            own = json.loads(p.read_text(encoding="utf-8"))
            p.with_name("scale_tile_estimate.json").write_text(json.dumps(own, indent=1), encoding="utf-8")
            new = {**own, "r_granulo": r_out, "openrec_radius": int(max(8, round(1.5 * r))),
                   "scale_rule": "per_scan", "scale_source": str(ROOT / "scale_central.json"),
                   "tile_estimate": {k: own.get(k) for k in ("r0", "r_granulo", "openrec_radius",
                                                              "specimen_fraction")}}
            p.write_text(json.dumps(new, indent=1), encoding="utf-8")
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(f"\n[fullvolume] scale rule per_scan: r_granulo {own.get('r_granulo')} -> {r_out}, "
                        f"openrec_radius {own.get('openrec_radius')} -> {new['openrec_radius']}\n")
        return rc, peak
    wrapped._gated = getattr(orig, "_gated", False)
    pl.run_tracked = wrapped
    log_note(f"scale rule per_scan installed: r_granulo {r_out}, openrec_radius {int(max(8, round(1.5 * r)))}")


def seen_stage_peaks(out: Path) -> dict:
    seen = {}
    for p in (out / "tiles").glob("*/tile.json"):
        r = json.loads(p.read_text(encoding="utf-8"))
        for s in r.get("steps") or []:
            k = Path(s["cmd0"]).name if s["cmd0"].endswith(".py") else {"scale": "batos_v3.scale",
                                                                           "flatten": "batos_v2.flatten_fixed"}.get(s["step"], s["cmd0"])
            seen[k] = max(seen.get(k, 0.0), float(s["peak_memory_mb"]))
    return seen


def write_manifest(name: str, tiles: list[Tile]) -> pd.DataFrame:
    rows = []
    for t in tiles:
        p = run_dir(name) / "tiles" / t.tile_id / "tile.json"
        if not p.exists():
            continue
        r = json.loads(p.read_text(encoding="utf-8"))
        rows.append({"tile_id": t.tile_id, "iz": t.iz, "iy": t.iy, "ix": t.ix,
                     **{f"core_{a}{e}": t.core[k][j] for k, a in enumerate("zyx") for j, e in enumerate("01")},
                     **{f"pad_{a}{e}": t.pad[k][j] for k, a in enumerate("zyx") for j, e in enumerate("01")},
                     **{k: r.get(k) for k in ("status", "r_granulo", "s", "size_filter_own", "openrec_radius",
                                              "openrec_impl", "stage03_thr", "specimen_fraction", "gate_triggered",
                                              "mv", "qv", "tos_markers", "BA-TOS_sel_instances", "BA-TOS_sel_owned",
                                              "BA-TOS_sel_owned_crossing_core", "v3_seconds", "post_seconds",
                                              "seconds_total", "peak_memory_mb", "git", "finished", "error")}})
    df = pd.DataFrame(rows)
    df.to_csv(run_dir(name) / "tiles_manifest.csv", index=False)
    return df


def cmd_run(a) -> int:
    global CORE
    CORE = int(a.core)
    cfg_method = check_freeze()
    sc = load_scale()
    halo = int(sc["halo"])
    shape = vol_shape()
    grid = make_grid(shape, CORE, halo)
    out = run_dir(a.name)
    (out / "tiles").mkdir(parents=True, exist_ok=True)
    sel_info = None
    if a.block == "specimen":
        tiles, sel_info, spec = specimen_selection(GRAY, grid)
        tifffile.imwrite(out / "specimen_mask_step8.tif", spec.astype(np.uint8), compression="zlib")
        fr = sel_info["padded_specimen_fraction"]
        # queue: specimen-filled tiles (padded coarse fraction >= 0.5) first in grid order, edge tiles last
        # by decreasing fraction (decided 2026-11-10)
        tiles = sorted(tiles, key=lambda t: (fr[t.tile_id] < 0.5, -fr[t.tile_id] if fr[t.tile_id] < 0.5 else 0,
                                             t.iz, t.iy, t.ix))
        sel_info["n_filled_ge050"] = int(sum(v >= 0.5 for v in fr.values()))
        log(f"specimen selection (GT-free): {sel_info}")
    else:
        tiles = select_tiles(grid, shape, CORE, a.block, a.tiles)
    rc_path = out / "run_config.json"
    run_cfg = {"name": a.name, "block": a.block, "halo": halo, "core": CORE, "volume_shape": list(shape),
               "tile_selection": sel_info, "scale_rule": a.scale_rule,
               "r_granulo_central": sc["r_granulo"], "method_tag": cfg_method["git_tag"],
               "method_code_commit": cfg_method["git_commit_code"], "fullvolume_config": CFG,
               "tiles": [{"tile_id": t.tile_id, "iz": t.iz, "iy": t.iy, "ix": t.ix, "core": t.core, "pad": t.pad}
                         for t in tiles]}
    if rc_path.exists():
        old = json.loads(rc_path.read_text(encoding="utf-8"))
        if ([t["tile_id"] for t in old["tiles"]] != [t.tile_id for t in tiles] or old["halo"] != halo
                or int(old["core"]) != CORE):
            raise RuntimeError(f"{rc_path} exists with a different tile set / halo; use another --name")
    else:
        rc_path.write_text(json.dumps(run_cfg, indent=1, default=list), encoding="utf-8")
    if a.dry_run:
        log(f"dry run: {len(tiles)} tiles, core {CORE}, halo {halo}, padded <= {CORE + 2 * halo}^3; "
            f"first {[t.tile_id for t in tiles[:3]]} last {[t.tile_id for t in tiles[-3:]]}")
        return 0
    sessions = out / "sessions.jsonl"
    sess = {"start": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
            "command": " ".join([Path(sys.executable).name, *sys.argv]), "git": git_commit(),
            "python": sys.version.split()[0], "platform": platform.platform(),
            "workers": a.workers, "mem_budget_gb": a.mem_budget_gb}
    todo = [t for t in tiles if not (out / "tiles" / t.tile_id / "tile.json").exists()]
    log(f"{a.name}: {len(tiles)} tiles ({len(tiles) - len(todo)} done), halo {halo}, padded <= "
        f"{CORE + 2 * halo}^3, workers <= {a.workers}, budget {a.mem_budget_gb} GB; frozen v3 OK "
        f"({cfg_method['git_tag']})")
    lock = threading.Lock()
    running: dict[str, threading.Thread] = {}
    errors: list[str] = []
    seen = seen_stage_peaks(out)
    for prev in a.seed_peaks_from or []:  # measured stage peaks of an earlier run (larger tiles = safe side)
        for k, v in seen_stage_peaks(run_dir(prev)).items():
            seen[k] = max(seen.get(k, 0.0), v)
    gate = StageGate(a.mem_budget_gb * 1024 - a.driver_overhead_mb, a.initial_peak_mb, seen)
    gate.install()
    if a.scale_rule == "per_scan":
        install_scale_rule(sc, log)
    log(f"stage gate: budget {gate.budget:.0f} MB for stage subprocesses; known stage peaks "
        f"{ {k: round(v) for k, v in gate.peak.items()} }")
    sess["stage_gate_budget_mb"] = gate.budget
    sampler = Sampler(out / "resources.csv", lambda: len(running))
    sampler.start()
    t_start = time.time()

    def work(t: Tile):
        try:
            process_tile(t, GRAY, out, halo)
        except Exception as e:  # driver bug, not a method failure: stop scheduling, report loudly
            errors.append(f"{t.tile_id}: {e!r}")
            log(f"ERROR {t.tile_id}: {e!r}")
        finally:
            with lock:
                running.pop(t.tile_id, None)

    queue = list(todo)
    while queue and not errors:
        with lock:
            n = len(running)
        if n < a.workers:
            t = queue.pop(0)
            th = threading.Thread(target=work, args=(t,), daemon=True)
            with lock:
                running[t.tile_id] = th
            log(f"start {t.tile_id} ({len(tiles) - len(queue)}/{len(tiles)}; running {n + 1}/{a.workers})")
            th.start()
        else:
            time.sleep(2)
    while running:
        time.sleep(2)
    sampler.stop_ev.set()
    sampler.join()
    sess |= {"end": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
             "wall_seconds": round(time.time() - t_start, 1), "tiles_run": len(todo) - len(queue),
             "tree_rss_peak_mb": round(sampler.peak_mb), "errors": errors}
    with open(sessions, "a", encoding="utf-8") as f:
        f.write(json.dumps(sess) + "\n")
    df = write_manifest(a.name, tiles)
    log(f"END {a.name}: wall {sess['wall_seconds']:.0f}s, tree RSS peak {sampler.peak_mb:.0f} MB, "
        f"status {df.status.value_counts().to_dict() if len(df) else {}}, errors {len(errors)}")
    return 1 if errors else 0


# ------------------------------------------------------------------------------------------------ stitch
def cmd_stitch(a) -> int:
    from batos_fullvolume.stitch import stitch
    tiles, cfg = load_tiles(a.name)
    missing = [t.tile_id for t in tiles if not (run_dir(a.name) / "tiles" / t.tile_id / "tile.json").exists()]
    if missing:
        raise RuntimeError(f"{len(missing)} tiles not finished, e.g. {missing[:5]}")
    meta = stitch(tiles, run_dir(a.name) / "tiles", run_dir(a.name) / f"stitch_{a.route}{a.variant}", a.route,
                  a.match_iou, orphan_rescue=a.orphan_rescue, rescue_min_free=a.rescue_min_free)
    log(f"stitched{a.variant} {meta['n_instances']} instances ({meta['n_dropped']} dropped, rescued "
        f"{meta['n_rescued']}, proxy seeds {meta['n_rescued_seed_proxy']}) in {meta['seconds']:.0f}s "
        f"box {meta['box']}")
    return 0


# ------------------------------------------------------------------------------------------------ evaluate
def gt_stats() -> pd.DataFrame:
    from batos_fullvolume.stream_eval import gt_label_stats
    p = ROOT / "gt_stats_EFRGP01_00.csv"
    if p.exists():
        return pd.read_csv(p)
    log("GT label statistics over the full volume (one streaming pass) ...")
    t0 = time.time()
    df = gt_label_stats(tifffile.memmap(GT, mode="r"), slab=16)
    df.to_csv(p, index=False)
    log(f"GT stats: {len(df)} labels in {time.time() - t0:.0f}s")
    return df


def roi_inside(path: Path = ROI) -> np.ndarray:
    vals = [int(x) for x in path.read_text(encoding="utf-8").split() if x.strip().lstrip("-").isdigit()]
    if not vals:
        raise RuntimeError(f"no labels read from {path}")
    return np.unique(np.asarray(vals, np.int64))


def crosses_grid(z0, z1, y0, y1, x0, x1):
    return ((z0 // CORE) != ((z1 - 1) // CORE)) | ((y0 // CORE) != ((y1 - 1) // CORE)) | \
        ((x0 // CORE) != ((x1 - 1) // CORE))


def cmd_evaluate(a) -> int:
    from batos_fullvolume.stream_eval import accumulate, finalize, seed_values
    t0 = time.time()
    tiles, cfg = load_tiles(a.name)
    sdir = run_dir(a.name) / f"stitch_{a.route}{a.variant}"
    meta = json.loads((sdir / "labels_meta.json").read_text(encoding="utf-8"))
    edir = run_dir(a.name) / f"eval_{a.route}{a.variant}{'' if a.gt_set == 'inside' else '_' + a.gt_set}"
    edir.mkdir(exist_ok=True)
    gs = gt_stats()
    roi = roi_inside()
    roi_out = roi_inside(ROI_OUT)
    core_box = [(min(t.core[k][0] for t in tiles), max(t.core[k][1] for t in tiles)) for k in range(3)]
    sel_cells = {(t.iz, t.iy, t.ix) for t in tiles}
    in_box = pd.Series([all((i, j, k) in sel_cells for i in range(z0 // CORE, (z1 - 1) // CORE + 1)
                            for j in range(y0 // CORE, (y1 - 1) // CORE + 1)
                            for k in range(x0 // CORE, (x1 - 1) // CORE + 1))
                        for z0, z1, y0, y1, x0, x1 in zip(gs.z0, gs.z1, gs.y0, gs.y1, gs.x0, gs.x1)],
                       index=gs.index)  # GT bbox inside the union of the selected tile cores
    gs["in_inside_list"] = gs.gt_label.isin(roi)
    gs["in_outside_list"] = gs.gt_label.isin(roi_out)
    gs["in_roi"] = gs.in_inside_list if a.gt_set == "inside" else (gs.in_inside_list | gs.in_outside_list)
    gs["valid"] = gs.in_roi & (gs.voxels >= 1000) & in_box
    valid_ids = np.sort(gs.loc[gs.valid, "gt_label"].to_numpy(np.int64))
    counts = {"gt_labels_full_volume": int(len(gs)), "gt_max_label": int(gs.gt_label.max()),
              "roi_list_labels": int(len(roi)), "roi_labels_present": int(gs.in_roi.sum()),
              "roi_present_ge1000": int((gs.in_roi & (gs.voxels >= 1000)).sum()),
              "gt_set": a.gt_set, "labels_inside_list": int(gs.in_inside_list.sum()),
              "labels_outside_list": int(gs.in_outside_list.sum()),
              "labels_unlisted": int((~gs.in_inside_list & ~gs.in_outside_list).sum()),
              "gt_fully_in_selected_cores": int(in_box.sum()), "valid_gt": int(len(valid_ids)), "core_box": core_box}
    log(f"GT counts: {counts}")
    L = np.load(sdir / "labels.npy", mmap_mode="r")
    gt = tifffile.memmap(GT, mode="r")
    off = tuple(meta["offset"])
    t1 = time.time()
    acc = accumulate(L, gt, valid_ids, pred_offset=off, slab=a.slab)
    sd = np.load(sdir / "seeds.npz")
    sp, sg = seed_values(L, gt, sd["zyx"], pred_offset=off)
    log(f"accumulated in {time.time() - t1:.0f}s; {len(sp)} seed voxels")
    inst = pd.read_csv(sdir / "instances.csv")
    ptot = acc["pred_total"]
    inst["voxels_final"] = ptot[inst.gid.to_numpy()] if len(inst) else []
    filters = {"none": np.zeros(0, np.int64),
               "500": inst.loc[inst.voxels_final < 500, "gid"].to_numpy(np.int64),
               "own": inst.loc[inst.voxels_final < inst.own_filter.fillna(500), "gid"].to_numpy(np.int64)}
    rows, details = [], {}
    for flt, rm in filters.items():
        m, table, det = finalize(acc, sp, sg, remove=rm)
        rows.append({"run": a.name, "route": a.route, "size_filter": flt, "n_removed_by_filter": int(len(rm)),
                     "mode": "markers_dominant", **m})
        details[flt] = (table, det)
        log(f"filter {flt}: F1 .25/.50/.75 {m['f1_iou025']:.4f}/{m['f1_iou050']:.4f}/{m['f1_iou075']:.4f} "
            f"P {m['precision_iou050']:.4f} R {m['recall_iou050']:.4f} pred_eval {m['prediction_count_evaluable']} "
            f"valid GT {m['valid_gt_count']}")
    met = pd.DataFrame(rows)
    met.to_csv(edir / "metrics.csv", index=False)
    # ---- per-GT / per-prediction detail at IoU 0.50, filter none (primary)
    table, det = details["none"]
    mp, mg, miou = det["matches"][0.5]
    gv = gs[gs.valid].copy()
    gv["owner_tile"] = [f"t{int(z) // CORE}_{int(y) // CORE}_{int(x) // CORE}"
                        for z, y, x in zip(gv.cz, gv.cy, gv.cx)]
    gv["gt_crosses_core"] = crosses_grid(gv.z0, gv.z1, gv.y0, gv.y1, gv.x0, gv.x1)
    gmatch = dict(zip(mg.tolist(), zip(mp.tolist(), miou.tolist())))
    gv["matched_gid"] = [gmatch.get(g, (0, np.nan))[0] for g in gv.gt_label]
    gv["matched_iou"] = [gmatch.get(g, (0, np.nan))[1] for g in gv.gt_label]
    # best-overlap evaluated prediction of each valid GT (for the reconciliation split of unmatched GT too)
    pp, gg, iou = det["pp"], det["gg"], det["iou"]
    best = pd.DataFrame({"p": det["eval_labels"][pp], "g": acc["valid_ids"][gg], "iou": iou})
    best = best.sort_values(["g", "iou"], ascending=[True, False]).drop_duplicates("g").set_index("g")
    gv["best_gid"] = gv.gt_label.map(best["p"]).fillna(0).astype(np.int64)
    gv["best_iou"] = gv.gt_label.map(best["iou"])
    if "rescued" not in inst:
        inst["rescued"] = 0
    inst["reconciled"] = ((inst.matched_replaced > 0) | (inst.contested_won > 0) | (inst.contested_lost > 0)
                          | (inst.rescued > 0))
    ii = inst.set_index("gid")
    gv["best_reconciled"] = gv.best_gid.map(ii["reconciled"]).fillna(False).astype(bool)
    gv["best_crosses_core"] = gv.best_gid.map(ii["crosses_core"]).fillna(False).astype(bool)
    gv.to_csv(edir / "per_gt.csv", index=False)
    ev = pd.DataFrame({"gid": det["eval_labels"]})
    ev = ev.merge(inst, on="gid", how="left")
    ev["matched"] = ev.gid.isin(set(mp.tolist()))
    ev.to_csv(edir / "per_pred.csv", index=False)
    table.to_csv(edir / "prediction_ignore_table.csv", index=False)
    # ---- per tile F1@0.50 (GT by centroid tile, predictions by owner tile)
    pt = []
    for t in tiles:
        g = gv[gv.owner_tile == t.tile_id]
        p = ev[ev.tile_id == t.tile_id]
        tp, fn, fp = int((g.matched_gid > 0).sum()), int((g.matched_gid == 0).sum()), int((~p.matched).sum())
        pt.append({"tile_id": t.tile_id, "iz": t.iz, "iy": t.iy, "ix": t.ix, "valid_gt": len(g), "pred_eval": len(p),
                   "tp": tp, "fp": fp, "fn": fn,
                   "f1_iou050": 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else np.nan})
    pd.DataFrame(pt).to_csv(edir / "per_tile.csv", index=False)

    # ---- seam check
    def f1_group(gmask, pmask):
        tp = int((gv.matched_gid[gmask] > 0).sum())
        fn = int((gv.matched_gid[gmask] == 0).sum())
        fp = int((~ev.matched[pmask]).sum())
        return {"valid_gt": int(gmask.sum()), "pred_eval": int(pmask.sum()), "tp": tp, "fp": fp, "fn": fn,
                "recall": tp / (tp + fn) if tp + fn else np.nan,
                "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else np.nan,
                "mean_matched_iou": float(gv.matched_iou[gmask & (gv.matched_gid > 0)].mean())}
    matched = gv.matched_gid > 0
    gv_rec = gv.matched_gid.map(ii["reconciled"]).fillna(False).astype(bool)
    gv_cross = gv.matched_gid.map(ii["crosses_core"]).fillna(False).astype(bool)
    seam = {
        "definitions": {
            "reconciled": "stitched instance that replaced / was compared to a matched instance of another tile "
                          "(IoU >= match_iou) or won / lost contested voxels to another tile's instance",
            "pred_crosses_core": "prediction extends beyond the core of the tile that owns it (carried by the halo)",
            "gt_crosses_core": "GT bbox spans >= 2 tile cores (seam GT); interior = inside one core",
            "group_f1": "TP / FN over the GT group, FP = unmatched evaluated predictions of the matching pred group"},
        "matched_gt": int(matched.sum()),
        "frac_matched_gt_pred_reconciled": float(gv_rec[matched].mean()) if matched.any() else np.nan,
        "frac_matched_gt_pred_crosses_core": float(gv_cross[matched].mean()) if matched.any() else np.nan,
        "frac_valid_gt_crosses_core": float(gv.gt_crosses_core.mean()) if len(gv) else np.nan,
        "n_instances": int(len(inst)), "n_reconciled_instances": int(inst.reconciled.sum()),
        "n_dropped_matched": int(len(pd.read_csv(sdir / "dropped.csv"))) if (sdir / "dropped.csv").stat().st_size > 2 else 0,
        "seam_gt_vs_interior": {
            "seam": f1_group(gv.gt_crosses_core.to_numpy(), ev.crosses_core.fillna(False).to_numpy(bool)),
            "interior": f1_group(~gv.gt_crosses_core.to_numpy(), ~ev.crosses_core.fillna(False).to_numpy(bool))},
        "reconciled_vs_not": {
            "reconciled": f1_group(gv.best_reconciled.to_numpy(), ev.reconciled.fillna(False).to_numpy(bool)),
            "not_reconciled": f1_group(~gv.best_reconciled.to_numpy(), ~ev.reconciled.fillna(False).to_numpy(bool))},
    }
    (edir / "seam.json").write_text(json.dumps(seam, indent=1, default=float), encoding="utf-8")
    # ---- resources
    man = pd.read_csv(run_dir(a.name) / "tiles_manifest.csv")
    sess = [json.loads(x) for x in (run_dir(a.name) / "sessions.jsonl").read_text(encoding="utf-8").splitlines() if x]
    res = {"tiles": len(man), "status": man.status.value_counts().to_dict(),
           "tile_seconds_sum": float(man.seconds_total.sum()), "tile_seconds_median": float(man.seconds_total.median()),
           "tile_seconds_max": float(man.seconds_total.max()), "wall_seconds_run": float(sum(s["wall_seconds"] for s in sess)),
           "tile_peak_mb_max": float(man.peak_memory_mb.max()), "tile_peak_mb_median": float(man.peak_memory_mb.median()),
           "tree_rss_peak_mb": float(max(s["tree_rss_peak_mb"] for s in sess)), "stitch_seconds": meta["seconds"],
           "eval_seconds": round(time.time() - t0, 1), "gate_triggered_tiles": int((man.gate_triggered == True).sum()),  # noqa: E712
           "r_granulo_tiles": man.r_granulo.describe().to_dict(), "r_granulo_central": cfg["r_granulo_central"],
           "halo": cfg["halo"]}
    out = {"run": a.name, "route": a.route, "variant": a.variant, "gt_set": a.gt_set, "core": CORE,
           "counts": counts, "resources": res, "stitch": meta,
           "metrics": rows, "seam": seam, "git": git_commit(),
           "timestamp": dt.datetime.now().astimezone().isoformat(timespec="seconds")}
    (edir / "metrics.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    log(f"seam: matched GT with reconciled pred {seam['frac_matched_gt_pred_reconciled']:.3f}; seam F1 "
        f"{seam['seam_gt_vs_interior']['seam']['f1']:.4f} vs interior {seam['seam_gt_vs_interior']['interior']['f1']:.4f}")
    log(f"evaluate done in {time.time() - t0:.0f}s -> {edir}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("scale")
    s.add_argument("--force", action="store_true")
    r = sub.add_parser("run")
    r.add_argument("--name", required=True)
    r.add_argument("--block", help="center<k> | all | specimen (GT-free specimen-touching tiles)")
    r.add_argument("--core", type=int, default=CORE, help="tile core edge (padded tile = core + 2 halo)")
    r.add_argument("--tiles", nargs="*")
    r.add_argument("--dry-run", action="store_true", help="select tiles, write run_config.json, run nothing")
    r.add_argument("--scale-rule", choices=["per_scan", "per_tile"],
                   default=CFG["tiling"].get("scale_rule", "per_scan"),
                   help="per_scan = central r_granulo for every tile (full-volume rule); per_tile = v3 as is")
    r.add_argument("--seed-peaks-from", nargs="*", help="earlier run names whose measured stage peaks seed the gate")
    r.add_argument("--workers", type=int, default=3)
    r.add_argument("--mem-budget-gb", type=float, default=float(CFG["resources"]["mem_budget_gb"]))
    r.add_argument("--initial-peak-mb", type=float, default=6000.0,
                   help="assumed per-STAGE peak before that stage has been measured")
    r.add_argument("--driver-overhead-mb", type=float, default=800.0,
                   help="RSS reserved for this driver process (post-processing of a 240^3 tile)")
    st = sub.add_parser("stitch")
    st.add_argument("--name", required=True)
    st.add_argument("--route", default=CFG["stitch"]["route"])
    st.add_argument("--match-iou", type=float, default=float(CFG["stitch"]["match_iou"]))
    st.add_argument("--orphan-rescue", action="store_true")
    st.add_argument("--rescue-min-free", type=float, default=0.5)
    st.add_argument("--variant", default="", help="suffix of the stitch dir, e.g. _rescue")
    e = sub.add_parser("evaluate")
    e.add_argument("--name", required=True)
    e.add_argument("--route", default=CFG["stitch"]["route"])
    e.add_argument("--slab", type=int, default=8)
    e.add_argument("--variant", default="", help="suffix of the stitch/eval dirs, e.g. _rescue")
    e.add_argument("--gt-set", choices=["inside", "listed"], default="inside",
                   help="inside = labelsInside (spec); listed = labelsInside + labelsOutside")
    a = ap.parse_args(argv)
    if getattr(a, "workers", 1) > int(CFG["resources"]["max_workers"]):
        raise SystemExit(f"--workers > {CFG['resources']['max_workers']} not allowed")
    return {"scale": cmd_scale, "run": cmd_run, "stitch": cmd_stitch, "evaluate": cmd_evaluate}[a.cmd](a)


if __name__ == "__main__":
    raise SystemExit(main())
