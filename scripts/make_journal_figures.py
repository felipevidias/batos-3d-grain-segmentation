#!/usr/bin/env python3
"""Qualitative journal figures (read-only on existing outputs; matplotlib only; 300 dpi).

Writes C:\\BA_TOS_DATA\\figures\\fig<N>_*.png and copies them to paper/figs/journal/. One figure per call
(--figs 1 2 ...) so peak memory stays bounded. Nothing is re-run: every segmentation, marker volume and scale estimate
is read from the frozen run folders below. F1 values are recomputed with the frozen evaluator
(batos_v2.evaluate_v2, no-marker fraction-0.50 ignore rule = the Cellpose comparison mode, no size filter) and
cross-checked against the recorded evaluations where those exist (mismatch -> abort).
"""

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
import tifffile  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from mpl_toolkits.axes_grid1.inset_locator import inset_axes  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
MAIN = Path(r"C:\Users\USER\batos-3d-grain-segmentation")  # results/cellpose lives on journal-extension only
sys.path.insert(0, str(REPO / "src"))
from batos_v2.evaluate_v2 import evaluate  # noqa: E402

D = Path(r"C:\BA_TOS_DATA")
R = D / "runs"
OUT = D / "figures"
PAPER = REPO / "paper" / "figs" / "journal"
DPI = 300
SEL = Path("e23c/watershed/segmentations/batos_selective_distance_c1.tif")
TOSM = Path("e11/watershed/normalized_markers/tos_e11_min.tif")
EVAL_KW = dict(rule="fraction", ignore_overlap_threshold=0.50, morph_iou=None)
METHODS = ["frozen", "fixflat", "v3", "cpsam_v2"]
MLABEL = {"frozen": "BA-TOS v2 (frozen)", "fixflat": "BA-TOS v2-fixflat", "v3": "BA-TOS v3",
          "cpsam_v2": "Cellpose-SAM v2"}
SPEC_COL = {"EFR": "#1f77b4", "EFS": "#17becf", "EMR": "#d62728", "EMS": "#ff9896", "EOR": "#2ca02c", "EOS": "#98df8a"}
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "legend.fontsize": 6,
                     "font.family": "DejaVu Sans"})


# ---------------------------------------------------------------- data access

def is_dev(cid: str) -> bool:
    return "_crop" in cid


def gray_gt(cid: str):
    base = D / ("crops" if is_dev(cid) else "crops_test") / cid
    g, l = base / "gray.tif", base / "labels.tif"
    if not (g.exists() and l.exists()):
        raise SystemExit(f"missing crop data {base}")
    return tifffile.imread(g), tifffile.imread(l)


def run_dir(method: str, cid: str) -> Path:
    if method == "cpsam_v2":
        return R / "cellpose" / "cpsam_v2" / cid
    if is_dev(cid):
        return {"frozen": R / "batos_paper_v2_frozen", "fixflat": R / "batos_paper_v2_fixflat" / "A",
                "v3": R / "v3_scale_test" / "runs" / "v3" / "x1"}[method] / cid
    return R / "test_protocol" / method / cid


def prediction(method: str, cid: str, shape) -> tuple[np.ndarray, str]:
    d = run_dir(method, cid)
    if (d / "FAILED.json").exists():
        return np.zeros(shape, np.uint32), "failed"
    p = d / ("labels.tif" if method == "cpsam_v2" else SEL)
    if not p.exists():
        raise SystemExit(f"missing prediction {p}")
    a = tifffile.imread(p)
    if a.shape != tuple(shape):
        raise SystemExit(f"shape mismatch {p}: {a.shape} vs {shape}")
    return a, "ok"


def recorded_f1(method: str, cid: str) -> float | None:
    """F1@0.50, BA-TOS_sel, no-marker fraction 0.50, no filter, from the evaluation of record (None if not stored)."""
    if method == "cpsam_v2":
        p = R / "cellpose" / "cpsam_v2" / "eval" / f"{cid}.csv"
        q = dict(size_filter="none")
    elif is_dev(cid) and method in ("v3", "fixflat"):
        p = R / "v3_scale_test" / "eval" / f"{cid}.csv"
        q = dict(method=method, route="BA-TOS_sel", size_filter="none", factor=1.0)
    elif not is_dev(cid):
        p = R / "test_protocol" / ("eval_v3" if method == "v3" else "eval") / f"{cid}.csv"
        q = dict(method=method, route="BA-TOS_sel", size_filter="none")
    else:
        return None
    if not p.exists():
        return None
    d = pd.read_csv(p, dtype={"size_filter": str})
    d = d[d["mode"] == "no_markers_fraction050"]
    for k, v in q.items():
        d = d[d[k] == v]
    return float(d.f1_iou050.iloc[0]) if len(d) == 1 else None


def f1_050(pred, gt, method, cid) -> float:
    m, _, _ = evaluate(pred, gt, None, **EVAL_KW)
    f = float(m["f1_iou050"])
    rec = recorded_f1(method, cid)
    if rec is not None and abs(rec - f) > 1e-6:
        raise SystemExit(f"F1 mismatch {method}/{cid}: recomputed {f:.6f} vs recorded {rec:.6f}")
    print(f"  {cid} {method}: F1@0.50 {f:.3f} (recorded {'n/a' if rec is None else f'{rec:.3f}'})", flush=True)
    return f


def gt_valid(gt: np.ndarray, min_vox: int = 1000) -> np.ndarray:
    """Boolean LUT over GT labels: >= min_vox voxels and not touching a crop face (evaluation protocol)."""
    n = int(gt.max()) + 1
    size = np.bincount(gt.ravel(), minlength=n)
    face = np.zeros(n, bool)
    for ax in range(3):
        for idx in (0, -1):
            face[np.unique(np.take(gt, idx, axis=ax))] = True
    ok = (size >= min_vox) & ~face
    ok[0] = False
    return ok


def match(pred: np.ndarray, gt: np.ndarray, thr: float = 0.5) -> dict[int, int]:
    """pred label -> GT label for 3D IoU > thr (unique by construction), plus pred -> dominant GT (key negative)."""
    p, g = pred.ravel().astype(np.int64), gt.ravel().astype(np.int64)
    sp, sg = np.bincount(p), np.bincount(g)
    both = (p > 0) & (g > 0)
    code = p[both] * (int(g.max()) + 1) + g[both]
    u, inter = np.unique(code, return_counts=True)
    pu, gu = u // (int(g.max()) + 1), u % (int(g.max()) + 1)
    iou = inter / (sp[pu] + sg[gu] - inter)
    out = {int(a): int(b) for a, b, i in zip(pu, gu, iou) if i > thr}
    order = np.argsort(-inter)  # dominant GT per pred (for dimming predictions inside excluded GT)
    dom: dict[int, int] = {}
    for a, b in zip(pu[order], gu[order]):
        dom.setdefault(int(a), int(b))
    return out | {-k: v for k, v in dom.items()}


def palette(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    h = rng.random(n)
    s = 0.55 + 0.45 * rng.random(n)
    v = 0.75 + 0.25 * rng.random(n)
    return matplotlib.colors.hsv_to_rgb(np.stack([h, s, v], 1))


def norm_gray(img: np.ndarray, ref: np.ndarray | None = None) -> np.ndarray:
    lo, hi = np.percentile(img if ref is None else ref, (1, 99.5))
    return np.clip((img.astype(np.float32) - lo) / max(hi - lo, 1e-6), 0, 1)


def overlay(gray2: np.ndarray, lab2: np.ndarray, cols: dict[int, np.ndarray], alpha: float = 0.65) -> np.ndarray:
    rgb = np.repeat(gray2[..., None], 3, 2)
    if lab2.max() == 0:
        return rgb
    lut = np.zeros((int(lab2.max()) + 1, 3), np.float32)
    for k, c in cols.items():
        if 0 < k < len(lut):
            lut[k] = c
    m = lab2 > 0
    rgb[m] = (1 - alpha) * rgb[m] + alpha * lut[lab2[m]]
    from skimage.segmentation import find_boundaries
    rgb[find_boundaries(lab2, mode="inner") & m] *= 0.35
    return rgb


def save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PAPER.mkdir(parents=True, exist_ok=True)
    p = OUT / name
    fig.savefig(p, dpi=DPI, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    shutil.copy2(p, PAPER / name)
    print(f"wrote {p} -> {PAPER / name}", flush=True)


# ---------------------------------------------------------------- fig 1

FIG1_CROPS = ["EFRGP01_00_crop04", "EFRGP01_00_crop05", "EMRGP01_00_test01", "EORGP01_00_test01"]


def fig1() -> None:
    from skimage.segmentation import find_boundaries  # noqa: F401
    cols_n = 2 + len(METHODS)
    fig, axs = plt.subplots(len(FIG1_CROPS), cols_n, figsize=(7.2, 1.38 * len(FIG1_CROPS) + 0.3))
    meta = []
    for r, cid in enumerate(FIG1_CROPS):
        gray, gt = gray_gt(cid)
        z = gt.shape[0] // 2
        g2 = norm_gray(gray[z], gray)
        ok = gt_valid(gt)
        gcol = palette(int(gt.max()) + 1, 7)
        dim = lambda c: 0.3 * c + 0.7 * np.array([0.55, 0.55, 0.55])  # noqa: E731
        gt_cols = {k: (gcol[k] if ok[k] else dim(gcol[k])) for k in np.unique(gt[z]) if k}
        panels = [("gray", g2, None, None), ("GT", overlay(g2, gt[z], gt_cols), None, None)]
        correct = []
        for mi, mth in enumerate(METHODS):
            pred, status = prediction(mth, cid, gt.shape)
            f1 = f1_050(pred, gt, mth, cid)
            mp = match(pred, gt)
            pcol = palette(int(pred.max()) + 1, 101 + mi)
            cols = {}
            for k in np.unique(pred[z]):
                if not k:
                    continue
                if k in mp:
                    cols[k] = gt_cols.get(mp[k], gcol[mp[k]] if ok[mp[k]] else dim(gcol[mp[k]]))
                elif -k in mp and not ok[mp[-k]]:
                    cols[k] = dim(pcol[k])
                else:
                    cols[k] = pcol[k]
            # per-pixel agreement with GT on the slice (for the inset window): matched label == GT label
            lut = np.zeros(int(pred.max()) + 1, np.int64)
            for k, v in mp.items():
                if k > 0:
                    lut[k] = v
            p2 = pred[z]
            correct.append(((lut[p2] == gt[z]) & (gt[z] > 0)) | ((p2 == 0) & (gt[z] == 0)))
            panels.append((mth, overlay(g2, p2, cols), f1, status))
            meta.append({"crop_id": cid, "method": mth, "status": status, "f1_050": f1, "slice_z": z})
            del pred
        # inset window: 48x48 where the methods disagree most (variance of per-pixel correctness), 10 px margin
        var = np.var(np.stack(correct).astype(np.float32), 0)
        W = 48
        sm = ndi.uniform_filter(var, W, mode="constant")
        sm[: W // 2 + 10], sm[-W // 2 - 10:], sm[:, : W // 2 + 10], sm[:, -W // 2 - 10:] = 0, 0, 0, 0
        cy, cx = np.unravel_index(np.argmax(sm), sm.shape)
        y0, x0 = cy - W // 2, cx - W // 2
        H = gt.shape[1]
        loc = ("lower " if cy < H / 2 else "upper ") + ("right" if cx < H / 2 else "left")  # opposite quadrant
        for c, (name, img, f1, status) in enumerate(panels):
            ax = axs[r, c]
            ax.imshow(img if img.ndim == 3 else img, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
            ax.set_xticks([]), ax.set_yticks([])
            if r == 0:
                ax.set_title({"gray": "grey", "GT": "ground truth"}.get(name, MLABEL.get(name, name)))
            if f1 is not None:
                ax.set_xlabel(f"F1 {f1:.3f}" + (" (run failed)" if status == "failed" else ""), labelpad=1.5)
            if c == 0:
                ax.set_ylabel(cid.replace("GP01_", " ").replace("_", " "), labelpad=2)
            ax.add_patch(Rectangle((x0 - 0.5, y0 - 0.5), W, W, fill=False, ec="yellow", lw=0.6))
            ins = inset_axes(ax, width="42%", height="42%", loc=loc, borderpad=0.15)
            ins.imshow(img[y0:y0 + W, x0:x0 + W], cmap="gray", vmin=0, vmax=1, interpolation="nearest")
            ins.set_xticks([]), ins.set_yticks([])
            for s in ins.spines.values():
                s.set_edgecolor("yellow"), s.set_linewidth(0.6)
        meta[-1]["inset_yx"] = [int(y0), int(x0), W]
        del gray, gt
    fig.subplots_adjust(wspace=0.03, hspace=0.12)
    save(fig, "fig1_qualitative_slices.png")
    pd.DataFrame(meta).to_csv(OUT / "fig1_f1.csv", index=False)


# ---------------------------------------------------------------- fig 2

def marker_centroids(path: Path) -> np.ndarray:
    m = tifffile.imread(path)
    lab = np.unique(m)
    lab = lab[lab > 0]
    if not len(lab):
        return np.zeros((0, 3))
    return np.array(ndi.center_of_mass(m > 0, m, lab))


def fig2() -> None:
    cid = "EFRGP01_00_crop05"
    gray, gt = gray_gt(cid)
    from skimage.filters import threshold_otsu
    from skimage.segmentation import find_boundaries
    sm = ndi.gaussian_filter(gray.astype(np.float32), 1.0)
    thr = threshold_otsu(sm)
    fg = sm > thr
    membrane = fg & (ndi.binary_dilation(gt > 0, iterations=1) == 0)
    lab, _ = ndi.label(membrane)
    if lab.max():
        membrane = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)  # the confining sheet (one component)
    dark, _ = ndi.label(~fg & ~ndi.binary_dilation(gt > 0, iterations=2))
    air = dark == (np.argmax(np.bincount(dark.ravel())[1:]) + 1)  # largest dark region = exterior air
    ok = gt_valid(gt)
    mk = {"frozen": marker_centroids(run_dir("frozen", cid) / TOSM), "v3": marker_centroids(run_dir("v3", cid) / TOSM)}
    print(f"  ToS-MIN markers: frozen {len(mk['frozen'])}, v3 {len(mk['v3'])}; air {air.mean():.3f}, "
          f"membrane {membrane.mean():.3f}", flush=True)
    c = np.array(gt.shape) // 2
    views = [("xy", 0, "x", "y"), ("xz", 1, "x", "z"), ("yz", 2, "y", "z")]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
    for ax, (name, axis, xl, yl) in zip(axs, views):
        sl = [slice(None)] * 3
        sl[axis] = c[axis]
        sl = tuple(sl)
        g2 = norm_gray(gray[sl], gray)
        rgb = np.repeat(g2[..., None], 3, 2)
        rgb[air[sl]] = 0.6 * rgb[air[sl]] + 0.4 * np.array([0.25, 0.45, 0.95])
        rgb[membrane[sl]] = 0.4 * rgb[membrane[sl]] + 0.6 * np.array([0.95, 0.3, 0.85])
        g = gt[sl]
        b = find_boundaries(g, mode="inner") & (g > 0)
        okl = ok[g]
        rgb[b & okl] = (0.1, 0.85, 0.1)
        rgb[b & ~okl] = (1.0, 0.6, 0.0)
        ax.imshow(rgb, interpolation="nearest")
        other = [i for i in range(3) if i != axis]
        cnt = {}
        for key, style in (("frozen", dict(c="red", marker="o", s=7, lw=0)),
                           ("v3", dict(facecolors="none", edgecolors="cyan", marker="o", s=14, lw=0.6))):
            p = mk[key]
            near = np.abs(p[:, axis] - c[axis]) <= 3 if len(p) else np.zeros(0, bool)
            ax.scatter(p[near, other[1]], p[near, other[0]], **style)
            cnt[key] = int(near.sum())
        ax.set_title(f"{name} slice ({'zyx'[axis]} = {c[axis]})\nmarkers within ±3: v2 {cnt['frozen']}, v3 {cnt['v3']}")
        ax.set_xlabel(xl), ax.set_ylabel(yl)
        ax.set_xticks([0, 100, 199]), ax.set_yticks([0, 100, 199])
    handles = [Patch(color=(0.1, 0.85, 0.1), label="valid GT boundary"),
               Patch(color=(1.0, 0.6, 0.0), label="excluded GT (face / < 1000 vox)"),
               Patch(color=(0.95, 0.3, 0.85), label="membrane (bright, unlabelled)"),
               Patch(color=(0.25, 0.45, 0.95), label="exterior air"),
               Line2D([], [], ls="", marker="o", ms=3, mfc="red", mec="red",
                      label=f"{MLABEL['frozen']} ToS-MIN markers ({len(mk['frozen'])} in crop)"),
               Line2D([], [], ls="", marker="o", ms=4, mfc="none", mec="cyan",
                      label=f"{MLABEL['v3']} ToS-MIN markers ({len(mk['v3'])} in crop)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.1))
    fig.subplots_adjust(wspace=0.32, bottom=0.17)
    save(fig, "fig2_crop05_orthogonal_markers.png")
    json.dump({"crop_id": cid, "otsu_smoothed": float(thr), "air_fraction": float(air.mean()),
               "membrane_fraction": float(membrane.mean()), "tos_min_markers": {k: len(v) for k, v in mk.items()}},
              open(OUT / "fig2_info.json", "w"), indent=2)


# ---------------------------------------------------------------- fig 3

def flat_path(variant: str, cid: str) -> Path:
    if variant == "frozen":
        d = R / "batos_paper_v2_frozen" / cid / "batos_openrec_auto_minparam" / "01_preprocess"
    else:
        d = R / "batos_paper_v2_fixflat" / "A" / cid / "pipe" / "01_preprocess"
    ps = sorted(d.glob("*_flatten_lambda1p0_8u.tif"))
    if len(ps) != 1:
        raise SystemExit(f"expected one flatten tif in {d}, got {ps}")
    return ps[0]


def fig3() -> None:
    cid = "EFRGP01_00_crop05"
    _, gt = gray_gt(cid)
    fs = pd.read_csv(MAIN / "results" / "flatten_fix" / "summary.csv")
    bg_fix = float(fs[(fs.crop_id == cid)].bg_ref.iloc[0])
    vols = {v: tifffile.imread(flat_path(v, cid)) for v in ("frozen", "fixflat")}
    ref = tifffile.imread(flat_path("frozen", "EFRGP01_00_crop04"))
    z = gt.shape[0] // 2
    fig, axs = plt.subplots(2, 2, figsize=(5.2, 4.6), gridspec_kw={"height_ratios": [1.0, 0.62]})
    axs[1, 1].sharey(axs[1, 0])
    titles = {"frozen": "frozen flatten (bg_ref = 0.9)", "fixflat": f"fixed flatten (bg_ref = {bg_fix:.1f})"}
    inside = gt > 0
    for j, (v, vol) in enumerate(vols.items()):
        ax = axs[0, j]
        ax.imshow(vol[z], cmap="gray", vmin=0, vmax=255, interpolation="nearest")
        ax.set_title(titles[v]), ax.set_xticks([]), ax.set_yticks([])
        h = axs[1, j]
        bins = np.arange(257)
        h.hist(vol[inside].ravel(), bins=bins, color="#2ca02c", alpha=0.75, label="inside GT grains", density=True)
        h.hist(vol[~inside].ravel(), bins=bins, color="#7f7f7f", alpha=0.6, label="outside (pores, membrane, air)",
               density=True)
        h.hist(ref[tifffile.imread(D / "crops" / "EFRGP01_00_crop04" / "labels.tif") > 0].ravel(), bins=bins,
               histtype="step", color="k", lw=0.7, label="crop04 grains (frozen; reference)", density=True)
        h.set_yscale("log"), h.set_xlim(0, 255), h.set_xlabel("selector grey value (8 bit)")
        if j == 0:
            h.set_ylabel("density (log)")
        med = float(np.median(vol[inside]))
        h.axvline(med, color="#2ca02c", ls="--", lw=0.7)
        h.text(med + 4, 0.93, f"grain median {med:.0f}", fontsize=6, color="#2ca02c",
               transform=h.get_xaxis_transform(), va="top")
        if j == 1:
            h.tick_params(labelleft=False)
    hh, ll = axs[1, 0].get_legend_handles_labels()
    fig.legend(hh, ll, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.04), fontsize=6)
    fig.subplots_adjust(hspace=0.12, wspace=0.06, bottom=0.12)
    save(fig, "fig3_crop05_flatten_offset.png")


# ---------------------------------------------------------------- fig 4

def fig4() -> None:
    dev = pd.read_csv(MAIN / "results" / "v3_scale_test" / "crop_info.csv")
    dev = dev[dev.factor == 1.0]
    dev = pd.DataFrame({"crop_id": dev.crop_id, "specimen": "EFR", "split": "dev", "gt_r": dev.gt_inscribed_r_median,
                        "r_granulo": dev.v3_r_scale, "r_e13": dev.fixflat_r_scale})
    t = pd.read_csv(MAIN / "results" / "test_v3" / "scale.csv")
    t = pd.DataFrame({"crop_id": t.crop_id, "specimen": t.specimen, "split": "test", "gt_r": t.gt_inscribed_r_median,
                      "r_granulo": t.v3_r_granulo, "r_e13": t.r_est_e13})
    d = pd.concat([dev, t], ignore_index=True)
    if len(d) != 64 or d[["gt_r", "r_granulo"]].isna().any().any():
        raise SystemExit(f"fig4: expected 64 crops with GT radius and r_granulo, got {len(d)} / NaN present")
    miss = d.loc[d.r_e13.isna(), "crop_id"].tolist()  # e13 failed upstream (recorded in test_protocol/crop_info.csv)
    print(f"  e13 r_est missing (e13 failed): {miss}", flush=True)
    d.to_csv(OUT / "fig4_scale.csv", index=False)
    lo, hi = 0, float(np.ceil(d[["gt_r", "r_granulo", "r_e13"]].max().max() + 1))
    fig, axs = plt.subplots(1, 2, figsize=(6.2, 3.0), sharey=True)
    rng = np.random.default_rng(0)
    for ax, (col, name) in zip(axs, (("r_granulo", r"$r_\mathrm{granulo}$ (v3)"), ("r_e13", r"$\hat r$ e13 (v2)"))):
        ax.fill_between([lo, hi], [lo * 0.9, hi * 0.9], [lo * 1.1, hi * 1.1], color="0.9", lw=0, label="±10 %")
        ax.plot([lo, hi], [lo, hi], "k-", lw=0.6)
        for sp, g in d.groupby("specimen"):
            for split, mk in (("dev", "s"), ("test", "o")):
                q = g[g.split == split]
                if len(q):
                    jx = rng.uniform(-0.08, 0.08, len(q))  # display jitter only; values are discrete
                    ax.scatter(q.gt_r + jx, q[col], s=16, marker=mk, color=SPEC_COL[sp], edgecolors="k",
                               linewidths=0.3, alpha=0.85, label=f"{sp} {split} ({len(q)})")
        q = d.dropna(subset=[col])
        rel = q[col] / q.gt_r - 1
        txt = (f"n = {len(q)}\nmedian rel. error {np.median(rel) * 100:+.0f} %\n"
               f"within ±10 %: {int((rel.abs() <= 0.10).sum())}/{len(q)}")
        if len(q) < len(d):
            txt += f"\n(e13 failed: {', '.join(d.loc[d[col].isna(), 'crop_id'])})"
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top", fontsize=6.5)
        ax.set_xlim(lo + 4, hi), ax.set_ylim(lo + 4, hi), ax.set_aspect("equal")
        ax.set_xlabel("GT median inscribed radius (voxels)"), ax.set_title(name)
    axs[0].set_ylabel("estimated grain radius (voxels)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="center left", bbox_to_anchor=(0.99, 0.5), frameon=False)
    save(fig, "fig4_scale_estimates.png")


# ---------------------------------------------------------------- fig 5

def fig5() -> None:
    from skimage.measure import find_contours
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 3.3))
    stats = []
    for ax, cid in zip(axs, ["EFRGP01_00_crop04", "EORGP01_00_test01"]):
        gray, gt = gray_gt(cid)
        ok = gt_valid(gt)
        z = gt.shape[0] // 2
        preds = {m: prediction(m, cid, gt.shape)[0] for m in ("v3", "cpsam_v2")}
        maps = {m: {v: k for k, v in match(p, gt).items() if k > 0} for m, p in preds.items()}  # GT -> pred
        common = [g for g in maps["v3"] if g in maps["cpsam_v2"] and ok[g]]
        sg = np.bincount(gt.ravel())
        rel = {}
        for m, p in preds.items():
            sp = np.bincount(p.ravel())
            rel[m] = np.array([sp[maps[m][g]] / sg[g] - 1 for g in common])
            stats.append({"crop_id": cid, "method": m, "n_common_matched": len(common),
                          "median_rel_volume_err": float(np.median(rel[m]))})
        # zoom window: 64x64 on the central slice with the most commonly matched valid grains
        W = 64
        on = np.isin(gt[z], common)
        cnt = ndi.uniform_filter(on.astype(np.float32), W, mode="constant")
        cnt[: W // 2], cnt[-W // 2:], cnt[:, : W // 2], cnt[:, -W // 2:] = 0, 0, 0, 0
        cy, cx = np.unravel_index(np.argmax(cnt), cnt.shape)
        y0, x0 = cy - W // 2, cx - W // 2
        win = (slice(y0, y0 + W), slice(x0, x0 + W))
        ax.imshow(norm_gray(gray[z], gray)[win], cmap="gray", vmin=0, vmax=1, interpolation="nearest")
        labs = [g for g in np.unique(gt[z][win]) if g in set(common)]
        for name, col, lw, src in (("GT", "white", 1.1, None), ("v3", "cyan", 0.8, "v3"),
                                   ("cpsam_v2", "magenta", 0.8, "cpsam_v2")):
            for g in labs:
                m2 = (gt[z][win] == g) if src is None else (preds[src][z][win] == maps[src][g])
                for cc in find_contours(np.pad(m2, 1).astype(float), 0.5):
                    ax.plot(cc[:, 1] - 1, cc[:, 0] - 1, color=col, lw=lw, ls="-" if src else "--")
        ax.set_xticks([]), ax.set_yticks([])
        ax.set_title(f"{cid}\nz = {z}, window y {y0}–{y0 + W}, x {x0}–{x0 + W}", fontsize=6.5)
        ax.set_xlabel(f"median volume error on {len(common)} grains matched by both\n"
                      f"v3 {np.median(rel['v3']) * 100:+.1f} %, Cellpose-SAM v2 {np.median(rel['cpsam_v2']) * 100:+.1f} %",
                      fontsize=6)
        del gray, gt, preds
    from matplotlib.lines import Line2D
    fig.legend([Line2D([], [], color="white", ls="--", lw=1.1, path_effects=None),
                Line2D([], [], color="cyan", lw=0.8), Line2D([], [], color="magenta", lw=0.8)],
               ["GT (dashed)", MLABEL["v3"], MLABEL["cpsam_v2"]], loc="lower center", ncol=3, frameon=True,
               facecolor="0.3", labelcolor="white", bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(wspace=0.08, bottom=0.16)
    save(fig, "fig5_boundary_contours.png")
    pd.DataFrame(stats).to_csv(OUT / "fig5_volume_error.csv", index=False)


# ---------------------------------------------------------------- fig 6

SMOKE = R / "fullvolume" / "smoke_center3" / "stitch_BA-TOS_sel_rescue"


def fig6() -> None:
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from skimage.measure import marching_cubes
    meta = json.loads((SMOKE / "labels_meta.json").read_text())
    lab = np.load(SMOKE / "labels.npy", mmap_mode="r")
    if list(lab.shape) != meta["shape"]:
        raise SystemExit("smoke labels shape mismatch")
    h = 20  # halo: evaluated core = box shrunk by the 20-voxel halo (run_config.json)
    # Ownership = instance centroid inside a tile core and no tiles exist outside this 3x3x3 block, so grains centred
    # outside it are absent and the outermost shell is half empty (labelled fraction 0.10-0.30 on the 600^3 faces vs
    # 0.50-0.64 inside). Panel (a) therefore cuts SHELL voxels in from every face.
    SHELL = 30
    core = lab[h + SHELL:-(h + SHELL), h + SHELL:-(h + SHELL), h + SHELL:-(h + SHELL)]
    n = core.shape[0]
    col = palette(int(meta["max_gid"]) + 1, 11)
    col[0] = (0.08, 0.08, 0.08)
    fig = plt.figure(figsize=(7.2, 3.6))
    # (a) three faces of the cut block, every 2nd voxel; faces shaded 1.0 / 0.85 / 0.7 for depth
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    st = 2
    idx = np.arange(0, n, st)
    A, B = np.meshgrid(idx, idx, indexing="ij")
    faces = [(1.0, np.asarray(core[n - 1, ::st, ::st]), (lambda a, b: (b, a, np.full_like(a, n - 1)))),
             (0.85, np.asarray(core[::st, 0, ::st]), (lambda a, b: (b, np.zeros_like(a), a))),
             (0.7, np.asarray(core[::st, ::st, n - 1]), (lambda a, b: (np.full_like(a, n - 1), b, a)))]
    for shade, f, xyz in faces:
        X, Y, Z = xyz(A, B)
        ax.plot_surface(X, Y, Z, facecolors=col[f] * shade, rstride=1, cstride=1, linewidth=0, antialiased=False,
                        shade=False)
    ax.set_box_aspect((1, 1, 1)), ax.view_init(elev=24, azim=-58), ax.set_axis_off()
    ax.set_title(f"(a) stitched BA-TOS v3 block (600³ core, 3×3×3 tiles),\nfaces cut {SHELL} voxels inside "
                 f"({n}³ shown)", fontsize=6.5)
    # (b) surfaces of the grains whose centroid lies in a central S^3 sub-block of the core (Lambert shading)
    ax = fig.add_subplot(1, 2, 2, projection="3d")
    S, M = 100, 40
    c0 = n // 2 - S // 2
    sub = np.asarray(core[c0 - M:c0 + S + M, c0 - M:c0 + S + M, c0 - M:c0 + S + M])
    ids = np.unique(sub)
    ids = ids[ids > 0]
    com = np.array(ndi.center_of_mass(np.ones_like(sub, bool), sub, ids)) if len(ids) else np.zeros((0, 3))
    keep = ids[((com >= M) & (com < M + S)).all(1)]
    objs = ndi.find_objects(sub)
    light = np.array([0.4, -0.5, 0.75])
    light /= np.linalg.norm(light)
    tris, fc, n_cut = [], [], 0
    for g in keep:
        sl = objs[g - 1]
        if any(s.start == 0 or s.stop == sub.shape[k] for k, s in enumerate(sl)):
            n_cut += 1  # touches the extraction margin (grain larger than M): skip rather than render clipped
            continue
        m = np.pad(sub[sl] == g, 1)
        if m.sum() < 500:
            continue
        v, f, _, _ = marching_cubes(m.astype(np.float32), 0.5, step_size=2)
        v = (v + np.array([s.start - 1 for s in sl]))[:, ::-1]  # zyx -> xyz
        t = v[f]
        nrm = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
        nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-9
        lam = 0.35 + 0.65 * np.abs(nrm @ light)
        tris.append(t)
        fc.append(col[g][None] * lam[:, None])
    pc = Poly3DCollection(np.concatenate(tris), facecolors=np.clip(np.concatenate(fc), 0, 1), linewidths=0)
    ax.add_collection3d(pc)
    lim = (M - 15, M + S + 15)
    ax.set_xlim(*lim), ax.set_ylim(*lim), ax.set_zlim(*lim)
    ax.set_box_aspect((1, 1, 1)), ax.view_init(elev=24, azim=-58), ax.set_axis_off()
    ax.set_title(f"(b) {len(tris)} grains with centroid in a central {S}³ sub-block\n(marching cubes)", fontsize=6.5)
    print(f"  fig6: {len(keep)} grains centred in sub-block, {n_cut} skipped (touch margin), {len(tris)} rendered")
    fig.subplots_adjust(wspace=0.0, left=0, right=1, top=0.93, bottom=0.0)
    save(fig, "fig6_smoke_block_3d.png")
    json.dump({"source": str(SMOKE), "core": 600, "shell_cut": SHELL, "shown": n, "n_instances": meta["n_instances"],
               "subblock": S, "subblock_origin_in_core": SHELL + c0, "grains_centred": int(len(keep)),
               "grains_skipped_margin": n_cut, "grains_rendered": len(tris)},
              open(OUT / "fig6_info.json", "w"), indent=2)


def main(argv=None) -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--figs", nargs="+", type=int, default=[1, 2, 3, 4, 5, 6])
    a = ap.parse_args(argv)
    for f in a.figs:
        print(f"fig {f}", flush=True)
        {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5, 6: fig6}[f]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
