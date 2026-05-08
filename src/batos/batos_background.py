#!/usr/bin/env python3
"""
filter_tos_seed_background_reject_v14.py

Patch v14 — correção direta para o erro observado:
as seeds ToSConOp estavam caindo no fundo/inter-grãos.

A ideia correta:
1. A ToSConOp continua gerando candidatos autoduais.
2. O volume flattened é usado para detectar regiões escuras.
3. Regiões escuras conectadas à BORDA DO VOLUME são classificadas como fundo externo.
4. A seed ToS só é aceita se NÃO estiver nesse fundo externo e se for compatível
   com um mínimo escuro interno cercado por uma vizinhança mais clara.

Isso é diferente dos filtros v8/v9/v13:
- não tenta adivinhar "grão inteiro" grosseiro;
- rejeita explicitamente fundo conectado à borda;
- mantém a ToSConOp como fonte dos candidatos.
"""

import argparse
from pathlib import Path
import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi


def read3d(path):
    arr = tiff.imread(str(path))
    arr = np.asarray(arr)
    if arr.ndim == 2:
        arr = arr[np.newaxis, :, :]
    if arr.ndim == 4 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    if arr.ndim != 3:
        raise RuntimeError(f"Esperava volume 3D, mas shape={arr.shape}")
    return arr


def robust_u8(vol, p_low=0.5, p_high=99.5):
    vol = np.asarray(vol)
    if vol.dtype == np.uint8:
        return vol.copy()
    v = vol.astype(np.float32)
    nz = v[v > 0]
    if nz.size == 0:
        nz = v.ravel()
    lo = np.percentile(nz, p_low)
    hi = np.percentile(nz, p_high)
    if hi <= lo:
        hi = lo + 1.0
    out = (v - lo) * 255.0 / (hi - lo)
    return np.clip(out, 0, 255).astype(np.uint8)


def otsu_u8(vol):
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


def choose_dark_threshold(flat_u8, dark_threshold, dark_percentile, otsu_offset, dark_max):
    if dark_threshold is not None and dark_threshold >= 0:
        return int(dark_threshold), otsu_u8(flat_u8), float(np.percentile(flat_u8[flat_u8 > 0], dark_percentile))

    nz = flat_u8[flat_u8 > 0]
    if nz.size == 0:
        nz = flat_u8.ravel()
    pval = float(np.percentile(nz, dark_percentile))
    otsu = otsu_u8(flat_u8)
    thr = min(otsu + otsu_offset, pval, dark_max)
    return int(np.clip(thr, 1, 255)), int(otsu), pval


def border_seed_from_mask(mask):
    seed = np.zeros_like(mask, dtype=bool)
    seed[0, :, :] |= mask[0, :, :]
    seed[-1, :, :] |= mask[-1, :, :]
    seed[:, 0, :] |= mask[:, 0, :]
    seed[:, -1, :] |= mask[:, -1, :]
    seed[:, :, 0] |= mask[:, :, 0]
    seed[:, :, -1] |= mask[:, :, -1]
    return seed


def compute_external_dark(flat_u8, dark_thr, connectivity=1):
    dark = flat_u8 <= int(dark_thr)
    border = border_seed_from_mask(dark)
    structure = ndi.generate_binary_structure(3, connectivity)
    external = ndi.binary_propagation(border, structure=structure, mask=dark)
    internal = dark & (~external)
    return dark, external, internal


def offsets_sphere(r):
    r = int(r)
    r2 = r * r
    offsets = []
    for dz in range(-r, r + 1):
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                d2 = dz*dz + dy*dy + dx*dx
                if d2 <= r2:
                    offsets.append((dz, dy, dx, d2))
    return offsets


def offsets_ring2d(inner, outer):
    i2 = int(inner) * int(inner)
    o2 = int(outer) * int(outer)
    r = int(outer)
    offsets = []
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            d2 = dx*dx + dy*dy
            if i2 < d2 <= o2:
                offsets.append((dy, dx, d2))
    return offsets


def sample_sphere(vol, z, y, x, offsets):
    D, H, W = vol.shape
    vals = []
    for dz, dy, dx, _ in offsets:
        zz, yy, xx = z + dz, y + dy, x + dx
        if 0 <= zz < D and 0 <= yy < H and 0 <= xx < W:
            vals.append(vol[zz, yy, xx])
    return np.asarray(vals) if vals else np.asarray([], dtype=vol.dtype)


def sample_ring2d(vol_or_mask, z, y, x, offsets):
    H, W = vol_or_mask.shape[1:]
    vals = []
    plane = vol_or_mask[z]
    for dy, dx, _ in offsets:
        yy, xx = y + dy, x + dx
        if 0 <= yy < H and 0 <= xx < W:
            vals.append(plane[yy, xx])
    return np.asarray(vals) if vals else np.asarray([], dtype=plane.dtype)


def find_near_internal(internal_dark, z, y, x, radius):
    if internal_dark[z, y, x]:
        return True, 0.0
    r = int(radius)
    D, H, W = internal_dark.shape
    best = None
    for dz in range(-r, r + 1):
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                zz, yy, xx = z + dz, y + dy, x + dx
                if not (0 <= zz < D and 0 <= yy < H and 0 <= xx < W):
                    continue
                if not internal_dark[zz, yy, xx]:
                    continue
                d2 = dz*dz + dy*dy + dx*dx
                if best is None or d2 < best:
                    best = d2
    if best is None:
        return False, np.inf
    return True, float(np.sqrt(best))


def to_rgb(gray):
    return np.repeat(gray[..., None], 3, axis=-1)


def overlay_mask(gray_u8, mask, color=(255, 70, 70), alpha=0.45):
    rgb = to_rgb(gray_u8).astype(np.float32)
    color = np.array(color, dtype=np.float32)
    m = mask > 0
    rgb[m] = (1.0 - alpha) * rgb[m] + alpha * color
    return np.clip(rgb, 0, 255).astype(np.uint8)


def overlay_external_internal(gray_u8, external_dark, internal_dark):
    rgb = to_rgb(gray_u8).astype(np.float32)
    e = external_dark > 0
    i = internal_dark > 0
    # fundo externo = vermelho
    rgb[e] = 0.55 * rgb[e] + 0.45 * np.array([255, 0, 0], dtype=np.float32)
    # mínimos internos = ciano
    rgb[i] = 0.55 * rgb[i] + 0.45 * np.array([0, 220, 255], dtype=np.float32)
    return np.clip(rgb, 0, 255).astype(np.uint8)


def overlay_centroids(gray_u8, centroids, size=1):
    rgb = to_rgb(gray_u8).astype(np.uint8)
    vals = np.unique(centroids)
    vals = vals[vals > 0]
    palette = [
        (255, 0, 0), (0, 255, 0), (0, 180, 255), (255, 255, 0),
        (255, 0, 255), (0, 255, 255), (255, 128, 0), (180, 0, 255),
        (255, 255, 255)
    ]
    D, H, W = centroids.shape
    for idx, v in enumerate(vals):
        coords = np.argwhere(centroids == v)
        if coords.size == 0:
            continue
        z, y, x = np.rint(coords.mean(axis=0)).astype(int)
        c = palette[idx % len(palette)]
        for dz in range(-size, size + 1):
            zz = z + dz
            if 0 <= zz < D:
                rgb[zz, y, x] = c
        for dy in range(-size, size + 1):
            yy = y + dy
            if 0 <= yy < H:
                rgb[z, yy, x] = c
        for dx in range(-size, size + 1):
            xx = x + dx
            if 0 <= xx < W:
                rgb[z, y, xx] = c
    return rgb


def nms_centroids(candidates, min_dist):
    if min_dist <= 0 or not candidates:
        return candidates, 0

    # prioriza candidatos com maior evidência interna e maior contraste
    candidates = sorted(candidates, key=lambda c: (-c["internal_overlap"], -c["ring_contrast"], c["area"]))
    accepted = []
    rejected = 0
    md2 = float(min_dist * min_dist)

    for c in candidates:
        p = np.array([c["z"], c["y"], c["x"]], dtype=np.float64)
        ok = True
        for a in accepted:
            q = np.array([a["z"], a["y"], a["x"]], dtype=np.float64)
            if np.sum((p - q) ** 2) < md2:
                ok = False
                break
        if ok:
            accepted.append(c)
        else:
            rejected += 1
    return accepted, rejected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("original")
    ap.add_argument("--flattened", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="bgreject")

    ap.add_argument("--dark-threshold", type=int, default=-1)
    ap.add_argument("--dark-percentile", type=float, default=38.0)
    ap.add_argument("--otsu-offset", type=int, default=-12)
    ap.add_argument("--dark-threshold-max", type=int, default=135)
    ap.add_argument("--dark-connectivity", type=int, default=1)

    ap.add_argument("--seed-area-min", type=int, default=1)
    ap.add_argument("--seed-area-max", type=int, default=2200)
    ap.add_argument("--max-external-overlap", type=float, default=0.01)
    ap.add_argument("--min-internal-overlap", type=float, default=0.02)
    ap.add_argument("--internal-search-radius", type=int, default=3)

    ap.add_argument("--center-radius", type=int, default=2)
    ap.add_argument("--ring-inner", type=int, default=4)
    ap.add_argument("--ring-outer", type=int, default=10)
    ap.add_argument("--center-dark-max", type=float, default=178.0)
    ap.add_argument("--min-ring-contrast", type=float, default=2.0)
    ap.add_argument("--max-ring-external-coverage", type=float, default=0.40)
    ap.add_argument("--nms-min-dist", type=float, default=4.0)

    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    original_u8 = robust_u8(read3d(args.original))
    flat_u8 = robust_u8(read3d(args.flattened))
    seed = read3d(args.seed)

    if flat_u8.shape != seed.shape:
        raise RuntimeError(f"shape mismatch: flat={flat_u8.shape}, seed={seed.shape}")

    dark_thr, otsu, pval = choose_dark_threshold(
        flat_u8,
        args.dark_threshold,
        args.dark_percentile,
        args.otsu_offset,
        args.dark_threshold_max,
    )

    dark_mask, external_dark, internal_dark = compute_external_dark(
        flat_u8,
        dark_thr,
        connectivity=args.dark_connectivity,
    )

    seed_lab, nseed = ndi.label(seed > 0, structure=ndi.generate_binary_structure(3, 1))
    slices = ndi.find_objects(seed_lab)

    center_sphere = offsets_sphere(args.center_radius)
    ring2d = offsets_ring2d(args.ring_inner, args.ring_outer)

    candidates = []
    reject = {
        "seed_area": 0,
        "external_background": 0,
        "not_internal_dark": 0,
        "center_not_dark": 0,
        "low_ring_contrast": 0,
        "ring_too_external": 0,
    }
    rows = []

    D, H, W = seed.shape

    for sid, slc in enumerate(slices, start=1):
        if slc is None:
            continue

        comp = seed_lab[slc] == sid
        area = int(comp.sum())

        if area < args.seed_area_min or area > args.seed_area_max:
            reject["seed_area"] += 1
            continue

        ext_overlap = float((comp & external_dark[slc]).sum()) / max(1, area)
        if ext_overlap > args.max_external_overlap:
            reject["external_background"] += 1
            rows.append((sid, 0, "external_background", area, ext_overlap, 0, 0, 0, 0, 0, 0))
            continue

        int_overlap = float((comp & internal_dark[slc]).sum()) / max(1, area)

        coords = np.argwhere(comp)
        z0, y0, x0 = slc[0].start, slc[1].start, slc[2].start
        cz, cy, cx = np.rint(coords.mean(axis=0) + np.array([z0, y0, x0])).astype(int)
        cz = int(np.clip(cz, 0, D - 1))
        cy = int(np.clip(cy, 0, H - 1))
        cx = int(np.clip(cx, 0, W - 1))

        near_internal, internal_dist = find_near_internal(internal_dark, cz, cy, cx, args.internal_search_radius)

        if int_overlap < args.min_internal_overlap and not near_internal:
            reject["not_internal_dark"] += 1
            rows.append((sid, 0, "not_internal_dark", area, ext_overlap, int_overlap, cz, cy, cx, internal_dist, 0))
            continue

        center_vals = sample_sphere(flat_u8, cz, cy, cx, center_sphere)
        ring_vals = sample_ring2d(flat_u8, cz, cy, cx, ring2d)
        ring_ext = sample_ring2d(external_dark.astype(np.uint8), cz, cy, cx, ring2d)

        center_mean = float(center_vals.mean()) if center_vals.size else 255.0
        ring_mean = float(ring_vals.mean()) if ring_vals.size else 0.0
        ring_contrast = ring_mean - center_mean
        ring_external_cov = float(ring_ext.mean()) if ring_ext.size else 1.0

        if center_mean > args.center_dark_max:
            reject["center_not_dark"] += 1
            rows.append((sid, 0, "center_not_dark", area, ext_overlap, int_overlap, cz, cy, cx, center_mean, ring_contrast))
            continue

        if ring_contrast < args.min_ring_contrast:
            reject["low_ring_contrast"] += 1
            rows.append((sid, 0, "low_ring_contrast", area, ext_overlap, int_overlap, cz, cy, cx, center_mean, ring_contrast))
            continue

        if ring_external_cov > args.max_ring_external_coverage:
            reject["ring_too_external"] += 1
            rows.append((sid, 0, "ring_too_external", area, ext_overlap, int_overlap, cz, cy, cx, ring_external_cov, ring_contrast))
            continue

        candidates.append({
            "sid": sid,
            "slice": slc,
            "component_mask": comp,
            "area": area,
            "z": cz,
            "y": cy,
            "x": cx,
            "external_overlap": ext_overlap,
            "internal_overlap": int_overlap,
            "center_mean": center_mean,
            "ring_mean": ring_mean,
            "ring_contrast": ring_contrast,
            "ring_external_cov": ring_external_cov,
        })

    kept_candidates, rejected_nms = nms_centroids(candidates, args.nms_min_dist)

    out_seed = np.zeros(seed.shape, dtype=np.uint8)
    centroids = np.zeros(seed.shape, dtype=np.uint16)

    for new_label, c in enumerate(kept_candidates, start=1):
        slc = c["slice"]
        comp = c["component_mask"]
        out_seed[slc][comp] = 255
        centroids[c["z"], c["y"], c["x"]] = np.uint16(new_label)
        rows.append((
            c["sid"], new_label, "keep", c["area"], c["external_overlap"], c["internal_overlap"],
            c["z"], c["y"], c["x"], c["center_mean"], c["ring_contrast"]
        ))

    # outputs
    tiff.imwrite(out_dir / f"{args.prefix}_dark_mask.tif", (dark_mask.astype(np.uint8) * 255))
    tiff.imwrite(out_dir / f"{args.prefix}_external_dark_background.tif", (external_dark.astype(np.uint8) * 255))
    tiff.imwrite(out_dir / f"{args.prefix}_internal_dark_candidates.tif", (internal_dark.astype(np.uint8) * 255))
    tiff.imwrite(out_dir / f"{args.prefix}_bgreject_seed_raw.tif", out_seed)
    tiff.imwrite(out_dir / f"{args.prefix}_bgreject_seed_centroids.tif", centroids)

    tiff.imwrite(out_dir / f"{args.prefix}_external_internal_overlay_rgb.tif",
                 overlay_external_internal(original_u8, external_dark, internal_dark))
    tiff.imwrite(out_dir / f"{args.prefix}_external_internal_overlay_on_flattened_rgb.tif",
                 overlay_external_internal(flat_u8, external_dark, internal_dark))
    tiff.imwrite(out_dir / f"{args.prefix}_bgreject_seed_overlay_rgb.tif",
                 overlay_mask(original_u8, out_seed, color=(255, 70, 70), alpha=0.45))
    tiff.imwrite(out_dir / f"{args.prefix}_bgreject_seed_overlay_on_flattened_rgb.tif",
                 overlay_mask(flat_u8, out_seed, color=(255, 70, 70), alpha=0.45))
    tiff.imwrite(out_dir / f"{args.prefix}_bgreject_centroids_overlay_rgb.tif",
                 overlay_centroids(original_u8, centroids, size=1))
    tiff.imwrite(out_dir / f"{args.prefix}_bgreject_centroids_overlay_on_flattened_rgb.tif",
                 overlay_centroids(flat_u8, centroids, size=1))

    with open(out_dir / f"{args.prefix}_bgreject_report.csv", "w") as f:
        f.write("seed_component,new_label,reason,seed_area,external_overlap,internal_overlap,z,y,x,metric1,metric2\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")

    print(f"[bgreject-v14] otsu={otsu}, percentile{args.dark_percentile}={pval:.3f}, dark_threshold={dark_thr}")
    print(f"[bgreject-v14] dark voxels     = {int(dark_mask.sum())}")
    print(f"[bgreject-v14] external dark   = {int(external_dark.sum())}")
    print(f"[bgreject-v14] internal dark   = {int(internal_dark.sum())}")
    print(f"[bgreject-v14] raw seed comps  = {nseed}")
    print(f"[bgreject-v14] candidates pre-NMS = {len(candidates)}")
    print(f"[bgreject-v14] rejected by NMS = {rejected_nms}")
    print(f"[bgreject-v14] kept final      = {len(kept_candidates)}")
    for k, v in reject.items():
        print(f"[bgreject-v14] rejected {k:22s} = {v}")
    print(f"[bgreject-v14] wrote = {out_dir}")
    print("[bgreject-v14] open first:")
    print(f"  gmic {out_dir}/{args.prefix}_external_internal_overlay_rgb.tif a z")
    print(f"  gmic {out_dir}/{args.prefix}_bgreject_seed_overlay_rgb.tif a z")
    print(f"  gmic {out_dir}/{args.prefix}_bgreject_centroids_overlay_rgb.tif a z")


if __name__ == "__main__":
    main()
