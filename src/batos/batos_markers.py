#!/usr/bin/env python3
"""
internal_dark_to_markers_batos.py

batos: quando a seed bruta da ToSConOp foi toda rejeitada como external_background,
não faz sentido continuar filtrando essa seed. O mapa ciano de "internal dark"
passa a ser a fonte de marcadores.

Entrada:
- original: crop original, só para overlay.
- flattened: volume flattened lambda=15.
- internal_dark: mapa binário de mínimos escuros internos gerado pelo v14
  (*_internal_dark_candidates.tif).

Saídas:
- *_batos_marker_components_raw.tif      componentes internos aceitos
- *_batos_marker_centroids.tif           centróides rotulados
- *_batos_marker_seed16_for_macke.tif    seed uint16 0/65535 para watershed.py da Macke
- overlays de QC
- CSV com métricas dos candidatos

Ideia:
internal dark candidates -> componentes 3D -> filtros geométricos/fotométricos -> NMS -> centróides.
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


def offsets_sphere(r):
    r = int(r)
    r2 = r * r
    out = []
    for dz in range(-r, r + 1):
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                d2 = dz*dz + dy*dy + dx*dx
                if d2 <= r2:
                    out.append((dz, dy, dx, d2))
    return out


def offsets_ring2d(inner, outer):
    i2 = int(inner) * int(inner)
    o2 = int(outer) * int(outer)
    r = int(outer)
    out = []
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            d2 = dx*dx + dy*dy
            if i2 < d2 <= o2:
                out.append((dy, dx, d2))
    return out


def sample_sphere(vol, z, y, x, offsets):
    D, H, W = vol.shape
    vals = []
    for dz, dy, dx, _ in offsets:
        zz, yy, xx = z + dz, y + dy, x + dx
        if 0 <= zz < D and 0 <= yy < H and 0 <= xx < W:
            vals.append(vol[zz, yy, xx])
    return np.asarray(vals) if vals else np.asarray([], dtype=vol.dtype)


def sample_ring2d(vol, z, y, x, offsets):
    H, W = vol.shape[1:]
    vals = []
    plane = vol[z]
    for dy, dx, _ in offsets:
        yy, xx = y + dy, x + dx
        if 0 <= yy < H and 0 <= xx < W:
            vals.append(plane[yy, xx])
    return np.asarray(vals) if vals else np.asarray([], dtype=vol.dtype)


def to_rgb(gray):
    return np.repeat(gray[..., None], 3, axis=-1)


def overlay_mask(gray_u8, mask, color=(255, 70, 70), alpha=0.45):
    rgb = to_rgb(gray_u8).astype(np.float32)
    color = np.asarray(color, dtype=np.float32)
    m = mask > 0
    rgb[m] = (1.0 - alpha) * rgb[m] + alpha * color
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


def nms(cands, min_dist, max_markers=0):
    if not cands:
        return [], 0

    # maior contraste primeiro; componentes muito grandes perdem prioridade
    cands = sorted(
        cands,
        key=lambda c: (-c["score"], c["area"], c["elongation"])
    )

    accepted = []
    rejected = 0
    md2 = float(min_dist * min_dist)

    for c in cands:
        p = np.asarray([c["z"], c["y"], c["x"]], dtype=np.float64)
        ok = True
        for a in accepted:
            q = np.asarray([a["z"], a["y"], a["x"]], dtype=np.float64)
            if float(np.sum((p - q) ** 2)) < md2:
                ok = False
                break
        if ok:
            accepted.append(c)
            if max_markers > 0 and len(accepted) >= max_markers:
                rejected += max(0, len(cands) - len(accepted))
                break
        else:
            rejected += 1

    return accepted, rejected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("original")
    ap.add_argument("--flattened", required=True)
    ap.add_argument("--internal-dark", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="batos")

    ap.add_argument("--connectivity", type=int, default=1)
    ap.add_argument("--area-min", type=int, default=3)
    ap.add_argument("--area-max", type=int, default=900)
    ap.add_argument("--max-dim", type=int, default=30)
    ap.add_argument("--max-elongation", type=float, default=8.0)
    ap.add_argument("--min-fill", type=float, default=0.005)

    ap.add_argument("--center-radius", type=int, default=2)
    ap.add_argument("--ring-inner", type=int, default=4)
    ap.add_argument("--ring-outer", type=int, default=10)
    ap.add_argument("--center-dark-max", type=float, default=205.0)
    ap.add_argument("--min-ring-contrast", type=float, default=-2.0)

    ap.add_argument("--nms-min-dist", type=float, default=7.0)
    ap.add_argument("--max-markers", type=int, default=420)

    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    original_u8 = robust_u8(read3d(args.original))
    flat_u8 = robust_u8(read3d(args.flattened))
    internal = read3d(args.internal_dark) > 0

    if flat_u8.shape != internal.shape:
        raise RuntimeError(f"shape mismatch: flat={flat_u8.shape}, internal={internal.shape}")

    structure = ndi.generate_binary_structure(3, args.connectivity)
    lab, n = ndi.label(internal, structure=structure)
    slices = ndi.find_objects(lab)

    center_sphere = offsets_sphere(args.center_radius)
    ring2d = offsets_ring2d(args.ring_inner, args.ring_outer)

    cands = []
    rows = []
    reject = {
        "area": 0,
        "shape": 0,
        "center": 0,
        "contrast": 0,
    }

    D, H, W = internal.shape

    for cid, slc in enumerate(slices, start=1):
        if slc is None:
            continue

        comp = lab[slc] == cid
        area = int(comp.sum())

        if area < args.area_min or area > args.area_max:
            reject["area"] += 1
            rows.append((cid, 0, "area", area, 0, 0, 0, 0, 0, 0, 0))
            continue

        dz = slc[0].stop - slc[0].start
        dy = slc[1].stop - slc[1].start
        dx = slc[2].stop - slc[2].start
        max_dim = max(dx, dy, dz)
        min_dim = max(1, min(dx, dy, dz))
        elong = max_dim / min_dim
        fill = area / max(1, dx * dy * dz)

        if max_dim > args.max_dim or elong > args.max_elongation or fill < args.min_fill:
            reject["shape"] += 1
            rows.append((cid, 0, "shape", area, dx, dy, dz, elong, fill, 0, 0))
            continue

        coords = np.argwhere(comp)
        z0, y0, x0 = slc[0].start, slc[1].start, slc[2].start
        cz, cy, cx = np.rint(coords.mean(axis=0) + np.asarray([z0, y0, x0])).astype(int)
        cz = int(np.clip(cz, 0, D - 1))
        cy = int(np.clip(cy, 0, H - 1))
        cx = int(np.clip(cx, 0, W - 1))

        center_vals = sample_sphere(flat_u8, cz, cy, cx, center_sphere)
        ring_vals = sample_ring2d(flat_u8, cz, cy, cx, ring2d)

        center_mean = float(center_vals.mean()) if center_vals.size else 255.0
        ring_mean = float(ring_vals.mean()) if ring_vals.size else 0.0
        ring_contrast = ring_mean - center_mean

        if center_mean > args.center_dark_max:
            reject["center"] += 1
            rows.append((cid, 0, "center", area, dx, dy, dz, elong, fill, center_mean, ring_contrast))
            continue

        if ring_contrast < args.min_ring_contrast:
            reject["contrast"] += 1
            rows.append((cid, 0, "contrast", area, dx, dy, dz, elong, fill, center_mean, ring_contrast))
            continue

        score = ring_contrast + max(0.0, 205.0 - center_mean) * 0.02

        cands.append({
            "component_id": cid,
            "slice": slc,
            "mask": comp,
            "area": area,
            "dx": dx,
            "dy": dy,
            "dz": dz,
            "elongation": elong,
            "fill": fill,
            "z": cz,
            "y": cy,
            "x": cx,
            "center_mean": center_mean,
            "ring_mean": ring_mean,
            "ring_contrast": ring_contrast,
            "score": score,
        })

    kept, rejected_nms = nms(cands, args.nms_min_dist, args.max_markers)

    seed_components = np.zeros(internal.shape, dtype=np.uint8)
    centroids = np.zeros(internal.shape, dtype=np.uint16)
    seed16 = np.zeros(internal.shape, dtype=np.uint16)

    for new_label, c in enumerate(kept, start=1):
        seed_components[c["slice"]][c["mask"]] = 255
        z, y, x = c["z"], c["y"], c["x"]
        centroids[z, y, x] = np.uint16(new_label)
        seed16[z, y, x] = np.uint16(65535)
        rows.append((
            c["component_id"], new_label, "keep", c["area"], c["dx"], c["dy"], c["dz"],
            c["elongation"], c["fill"], c["center_mean"], c["ring_contrast"]
        ))

    tiff.imwrite(out_dir / f"{args.prefix}_batos_marker_components_raw.tif", seed_components)
    tiff.imwrite(out_dir / f"{args.prefix}_batos_marker_centroids.tif", centroids)
    tiff.imwrite(out_dir / f"{args.prefix}_batos_marker_seed16_for_macke.tif", seed16)

    tiff.imwrite(out_dir / f"{args.prefix}_batos_internal_dark_overlay_rgb.tif",
                 overlay_mask(original_u8, internal, color=(0, 220, 255), alpha=0.40))
    tiff.imwrite(out_dir / f"{args.prefix}_batos_marker_components_overlay_rgb.tif",
                 overlay_mask(original_u8, seed_components, color=(255, 70, 70), alpha=0.45))
    tiff.imwrite(out_dir / f"{args.prefix}_batos_marker_components_overlay_on_flattened_rgb.tif",
                 overlay_mask(flat_u8, seed_components, color=(255, 70, 70), alpha=0.45))
    tiff.imwrite(out_dir / f"{args.prefix}_batos_marker_centroids_overlay_rgb.tif",
                 overlay_centroids(original_u8, centroids, size=1))
    tiff.imwrite(out_dir / f"{args.prefix}_batos_marker_centroids_overlay_on_flattened_rgb.tif",
                 overlay_centroids(flat_u8, centroids, size=1))

    with open(out_dir / f"{args.prefix}_batos_marker_report.csv", "w") as f:
        f.write("component_id,new_label,reason,area,dx,dy,dz,elongation,fill,center_mean,ring_contrast\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")

    print(f"[BA-TOS] internal dark components = {n}")
    print(f"[BA-TOS] candidates before NMS    = {len(cands)}")
    print(f"[BA-TOS] rejected by NMS          = {rejected_nms}")
    print(f"[BA-TOS] kept markers             = {len(kept)}")
    for k, v in reject.items():
        print(f"[BA-TOS] rejected {k:10s} = {v}")
    print(f"[BA-TOS] wrote = {out_dir}")
    print("[BA-TOS] open:")
    print(f"  gmic {out_dir}/{args.prefix}_batos_internal_dark_overlay_rgb.tif a z")
    print(f"  gmic {out_dir}/{args.prefix}_batos_marker_components_overlay_rgb.tif a z")
    print(f"  gmic {out_dir}/{args.prefix}_batos_marker_centroids_overlay_rgb.tif a z")


if __name__ == "__main__":
    main()
