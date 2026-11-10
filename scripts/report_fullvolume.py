#!/usr/bin/env python3
"""Phase 5g report: copies the evaluation tables of one full-volume run to results/fullvolume/<run>/, draws the
per-tile F1@0.50 heat map (one panel per z-slab, tile grid y x x) and upserts pooled rows into results/all_runs.csv.
Pooled-object metrics only (a full volume has no crop-macro mean; the per-tile F1 is a spatial diagnostic)."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
ROOT = Path(r"C:\BA_TOS_DATA\runs\fullvolume")
# single-hue sequential ramp (dataviz reference palette, blue 100 -> 700)
SEQ = LinearSegmentedColormap.from_list("blue_seq", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])


def upsert(path: Path, rows: list[dict], keys: list[str]) -> None:
    new = pd.DataFrame(rows)
    if path.exists():
        old = pd.read_csv(path)
        k_old = old[keys].astype(str).agg("|".join, axis=1)
        k_new = set(new[keys].astype(str).agg("|".join, axis=1))
        old = old[~k_old.isin(k_new)]
        new = pd.concat([old, new], ignore_index=True)[old.columns.union(new.columns, sort=False)]
    new.to_csv(path, index=False)


def heatmap(pt: pd.DataFrame, out: Path, title: str, core: int = 200) -> None:
    zs = sorted(pt.iz.unique())
    ys, xs = sorted(pt.iy.unique()), sorted(pt.ix.unique())
    n = len(zs)
    cols = min(n, 4)
    rows = int(np.ceil(n / cols))
    cell = 0.55 if len(xs) > 4 else 0.9
    fig, axes = plt.subplots(rows, cols, figsize=(cols * len(xs) * cell + 1.6, rows * len(ys) * cell + 1.2),
                             squeeze=False, constrained_layout=True)
    vmin = max(0.0, np.floor(np.nanmin(pt.f1_iou050.to_numpy(float)) * 10) / 10) if pt.f1_iou050.notna().any() else 0
    im = None
    for k, z in enumerate(zs):
        ax = axes[k // cols][k % cols]
        g = np.full((len(ys), len(xs)), np.nan)
        ngt = np.zeros_like(g)
        for r in pt[pt.iz == z].itertuples():
            g[ys.index(r.iy), xs.index(r.ix)] = r.f1_iou050
            ngt[ys.index(r.iy), xs.index(r.ix)] = r.valid_gt
        ax.set_facecolor("#f0efec")
        im = ax.imshow(g, cmap=SEQ, vmin=vmin, vmax=1.0, origin="upper")
        for i in range(len(ys)):
            for j in range(len(xs)):
                if np.isfinite(g[i, j]):
                    ax.text(j, i, f"{g[i, j]:.3f}\nn={int(ngt[i, j])}", ha="center", va="center", fontsize=7,
                            color="white" if g[i, j] > vmin + 0.6 * (1 - vmin) else "#1f1f1f")
                else:
                    ax.text(j, i, "no GT", ha="center", va="center", fontsize=6, color="#6b6b6b")
        ax.set_xticks(range(len(xs)), [f"x{v}" for v in xs], fontsize=7)
        ax.set_yticks(range(len(ys)), [f"y{v}" for v in ys], fontsize=7)
        ax.set_title(f"z-slab {z} (z {z * core}-{z * core + core})", fontsize=8)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
    for k in range(n, rows * cols):
        axes[k // cols][k % cols].axis("off")
    cb = fig.colorbar(im, ax=axes, shrink=0.8)
    cb.set_label("F1@0.50 (valid GT by centroid tile; FP by owner tile)", fontsize=8)
    cb.outline.set_visible(False)
    fig.suptitle(title, fontsize=9)
    fig.savefig(out, dpi=160)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--name", required=True)
    ap.add_argument("--route", default="BA-TOS_sel")
    ap.add_argument("--variant", default="", help="stitch/eval dir suffix, e.g. _rescue")
    ap.add_argument("--gt-set", choices=["inside", "listed"], default="inside")
    a = ap.parse_args(argv)
    gsfx = "" if a.gt_set == "inside" else "_" + a.gt_set
    edir = ROOT / a.name / f"eval_{a.route}{a.variant}{gsfx}"
    res = REPO / "results" / "fullvolume" / f"{a.name}{a.variant}{gsfx}"
    res.mkdir(parents=True, exist_ok=True)
    for f in ("metrics.csv", "metrics.json", "per_tile.csv", "seam.json"):
        shutil.copy2(edir / f, res / f)
    shutil.copy2(ROOT / a.name / "tiles_manifest.csv", res / "tiles_manifest.csv")
    js = json.loads((edir / "metrics.json").read_text(encoding="utf-8"))
    pt = pd.read_csv(edir / "per_tile.csv")
    core = int(js.get("core", 200))
    heatmap(pt, res / "per_tile_f1_heatmap.png", core=core, title=
            f"{a.name}{a.variant}{gsfx}: BA-TOS v3 tiled ({a.route}, seed-ignore, no size filter) - per-tile F1@0.50")
    man = pd.read_csv(ROOT / a.name / "tiles_manifest.csv")
    r = js["resources"]
    rows = []
    for m in js["metrics"]:
        rows.append({"experiment_id": f"fullvolume_{a.name}{a.variant}{gsfx}/EFRGP01_00", "method": a.route,
                     "model": "batos-paper-v3-scalefree (tiled)", "volume": "EFRGP01_00",
                     "crop": f"{a.name} ({r['tiles']} tiles)", "device": "cpu (AMD Ryzen 7 5700X)",
                     "parameters": f"size_filter={m['size_filter']};eval=stream_eval:markers_dominant;"
                                   f"valid={'roi_inside' if a.gt_set == 'inside' else 'roi_inside+outside'}+ge1000+in_selected_cores;"
                                   f"core={core};halo={r['halo']};"
                                   f"match_iou={js['stitch']['match_iou']};contested=core_authority;"
                                   f"orphan_rescue={js['stitch'].get('orphan_rescue', False)}",
                     "runtime_seconds": r["tile_seconds_sum"], "peak_memory_mb": r["tree_rss_peak_mb"],
                     "tp": m["tp_iou050"], "fp": m["fp_iou050"], "fn": m["fn_iou050"],
                     "precision": m["precision_iou050"], "recall": m["recall_iou050"], "f1_025": m["f1_iou025"],
                     "f1_050": m["f1_iou050"], "f1_075": m["f1_iou075"], "ap50": m["ap_iou050"],
                     "ap75": m["ap_iou075"], "mean_matched_iou": m["mean_matched_iou050"],
                     "n_pred": m["prediction_count_evaluable"], "n_gt": m["valid_gt_count"],
                     "git_commit": js["git"]["commit"], "run_json": str(edir / "metrics.json"),
                     "memory_run_json": str(ROOT / a.name / "resources.csv")})
    upsert(REPO / "results" / "all_runs.csv", rows, ["experiment_id", "method", "parameters"])
    print(f"report -> {res}; {len(rows)} rows upserted; tiles {len(man)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
