#!/usr/bin/env python3
"""Phase 1f: define the TEST crops by a fixed, GT-free rule and write data/metadata/crops_test.json.

Run BEFORE any test crop is cut. Reads ONLY the greyscale volumes (memory-mapped, every 4th voxel); never the GT,
never the ROI label lists. Six Pinzon et al. specimens (EFR, EFS, EMR, EMS, EOR, EOS), scans _00 and _01.
EFR = internal test (6 crops per scan, no overlap with the 12 development crops); the five others = external test
(4 crops per scan). 52 crops in total.

Rule (also stored verbatim in the JSON):
 1. Specimen mask per scan: greyscale subsampled by 4 on every axis; one Otsu threshold over all nonzero voxels;
    per z-slice: binary closing (3 iterations), fill holes, largest 2-D component. This is the specimen cross-section
    bounded by the confining membrane (the thin bright sheet between grains and air); beyond it is exterior air.
 2. Specimen z-range: slices whose equivalent radius sqrt(area/pi) is within +-3 % of the median, in BOTH scans of the
    specimen (the end caps / platens have a larger cross-section). Box z-extent [z, z+200) must lie inside it.
 3. Specimen axis: median over those slices of the cross-section centroid, per scan.
 4. Candidates: origins on a 25-voxel lattice (x, y, z), box fully inside the volume and the specimen z-range.
    A box may not overlap (any voxel) an already selected test box of the same specimen (the frame is shared by the
    two scans, so this also holds across scans), nor, for EFR, any of the 6 development origins (12 dev crops).
 5. f = specimen fraction of the box (mean of the subsampled mask over the box, in the scan the box is cut from).
    Class "interior": f >= 0.99. Class "boundary" (intersects the membrane, contains exterior air): 0.40 <= f <= 0.80
    (the EFR development crops span f = 0.50 (crop05) .. 0.79 (crop03) .. 1.00).
 6. Strata per specimen: z-tier (origin z in equal parts of the admissible range) x angular half of the box centre
    around the specimen axis (theta = atan2(dy, dx) in [0, pi) -> scan _00, [-pi, 0) -> scan _01) x class
    (interior / boundary). One crop per cell. EFR: 3 z-tiers (3 interior + 3 boundary per scan); external specimens:
    2 z-tiers (2 interior + 2 boundary per scan). Every scan has >= 1 crop that intersects the specimen boundary.
 7. Cells are filled in a fixed order (specimen EFR, EFS, EMR, EMS, EOR, EOS; scan half; tier; class); in each cell one
    admissible candidate is drawn uniformly at random with ONE numpy.random.default_rng(SEED) shared by all specimens.
    If a cell has no candidate the script fails (no relaxation).
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu

REPO = Path(__file__).resolve().parents[1]
DEV = REPO / "data" / "metadata" / "crops_sibgrapi.json"
OUT = REPO / "data" / "metadata" / "crops_test.json"
DATA = Path(r"C:\Users\USER\Downloads\PhDUpload\PhDUpload")
CROPS_ROOT = Path(r"C:\BA_TOS_DATA\crops_test")
SEED = 20261106
SUB = 4
SIZE = 200
LATTICE = 25
F_INTERIOR = 0.99
F_BOUNDARY = (0.40, 0.80)
RADIUS_TOL = 0.03
SPECIMENS = {"EFR": 3, "EFS": 2, "EMR": 2, "EMS": 2, "EOR": 2, "EOS": 2}  # z-tiers; crops per scan = 2 x tiers
CLASSES = ("interior", "boundary")


def specimen_mask(gray_path: str):
    vol = tifffile.memmap(gray_path, mode="r")
    sub = np.array(vol[::SUB, ::SUB, ::SUB])
    th = float(threshold_otsu(sub[sub > 0]))
    mask = np.zeros(sub.shape, bool)
    for z in range(sub.shape[0]):
        fill = ndi.binary_fill_holes(ndi.binary_closing(sub[z] > th, iterations=3))
        lab, n = ndi.label(fill)
        if n:
            sizes = np.bincount(lab.ravel())
            sizes[0] = 0
            mask[z] = lab == sizes.argmax()
    req = np.sqrt(mask.reshape(mask.shape[0], -1).sum(1) / np.pi) * SUB
    ok = np.abs(req / np.median(req) - 1) <= RADIUS_TOL
    return vol.shape, mask, th, req, ok


def geometry(gray_path: str) -> dict:
    shape, mask, th, req, ok = specimen_mask(gray_path)
    zs = np.nonzero(ok)[0]
    cy, cx = [], []
    for z in zs:
        yy, xx = np.nonzero(mask[z])
        cy.append(yy.mean() * SUB)
        cx.append(xx.mean() * SUB)
    sat = np.zeros(tuple(n + 1 for n in mask.shape), np.int64)
    sat[1:, 1:, 1:] = mask.cumsum(0).cumsum(1).cumsum(2)
    return {"shape": [int(v) for v in shape], "otsu": th, "z_ok": [int(zs.min()) * SUB, int(zs.max()) * SUB],
            "axis_xy": [float(np.median(cx)), float(np.median(cy))], "req_median": float(np.median(req)),
            "mask_sha256": hashlib.sha256(np.packbits(mask).tobytes()).hexdigest(), "_sat": sat}


def frac(g: dict, c: dict) -> float:
    """Specimen fraction of the box: mean of the subsampled mask over [z//4, (z+200)//4) x ... (summed-area table)."""
    s, n = g["_sat"], SIZE // SUB
    z0, y0, x0 = c["z"] // SUB, c["y"] // SUB, c["x"] // SUB
    z1, y1, x1 = z0 + n, y0 + n, x0 + n
    tot = (s[z1, y1, x1] - s[z0, y1, x1] - s[z1, y0, x1] - s[z1, y1, x0]
           + s[z0, y0, x1] + s[z0, y1, x0] + s[z1, y0, x0] - s[z0, y0, x0])
    return float(tot) / n ** 3


def overlaps(a, b) -> bool:
    return all(abs(a[k] - b[k]) < SIZE for k in ("x", "y", "z"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    if args.out.exists():
        raise FileExistsError(f"{args.out} exists; the test set is defined once")
    dev = json.loads(DEV.read_text(encoding="utf-8"))
    dev_boxes = [dict(x=x, y=y, z=z) for x, y, z in sorted({(c["x"], c["y"], c["z"]) for c in dev["crops"]})]
    rng = np.random.default_rng(args.seed)
    chosen, geo_out, sources, params_spec, dev_f = [], {}, {}, {}, {}
    for spec, n_tiers in SPECIMENS.items():
        scans = {0: f"{spec}GP01_00", 1: f"{spec}GP01_01"}
        geo = {}
        for scan in scans.values():
            gray = DATA / spec / "Greyscale" / f"{scan}.tif"
            geo[scan] = geometry(str(gray))
            g = geo[scan]
            sources[scan] = {"gray_path": str(gray), "shape_zyx": g["shape"], "dtype": "uint16",
                             "labels_archive": str(DATA / spec / "Label" / f"{spec}_labels.7z"),
                             "labels_member": f"{scan}.tif"}
            print(f"[geo] {scan}: otsu {g['otsu']:.0f}, R_eq {g['req_median']:.1f}, z_ok {g['z_ok']}, "
                  f"axis x,y {g['axis_xy'][0]:.1f},{g['axis_xy'][1]:.1f}", flush=True)
        shape = geo[scans[0]]["shape"]
        z_lo = max(g["z_ok"][0] for g in geo.values())
        z_hi = min(g["z_ok"][1] for g in geo.values())  # last admissible slice (subsampled grid, inclusive)
        z_min = int(np.ceil(z_lo / LATTICE) * LATTICE)
        z_max = int(np.floor((z_hi + SUB - SIZE) / LATTICE) * LATTICE)  # z + SIZE <= z_hi + SUB
        tiers = np.linspace(z_min, z_max + LATTICE, n_tiers + 1)
        params_spec[spec] = {"admissible_origin_z": [z_min, z_max], "z_tier_edges": tiers.tolist(),
                             "scan_of_angular_half": scans,
                             "excluded_boxes": "12 development crops (6 origins)" if spec == "EFR" else "none"}
        print(f"[geo] {spec}: admissible origin z {z_min}..{z_max}; tier edges {tiers.tolist()}", flush=True)
        excluded = dev_boxes if spec == "EFR" else []
        lattice = [dict(x=x, y=y, z=z) for z in range(z_min, z_max + 1, LATTICE)
                   for y in range(0, shape[1] - SIZE + 1, LATTICE) for x in range(0, shape[2] - SIZE + 1, LATTICE)]
        picked = []
        for half, scan in scans.items():
            ax, ay = geo[scan]["axis_xy"]
            for tier in range(n_tiers):
                for cls in CLASSES:
                    cands = []
                    for c in lattice:
                        if not (tiers[tier] <= c["z"] < tiers[tier + 1]):
                            continue
                        theta = np.arctan2(c["y"] + SIZE / 2 - ay, c["x"] + SIZE / 2 - ax)
                        if (0 if theta >= 0 else 1) != half:
                            continue
                        if any(overlaps(c, b) for b in excluded + picked):
                            continue
                        f = frac(geo[scan], c)
                        if (cls == "interior" and f >= F_INTERIOR) or \
                                (cls == "boundary" and F_BOUNDARY[0] <= f <= F_BOUNDARY[1]):
                            cands.append((c, f, theta))
                    if not cands:
                        raise RuntimeError(f"no admissible candidate for {scan} tier {tier} {cls}")
                    c, f, theta = cands[int(rng.integers(len(cands)))]
                    n = len(picked) + 1
                    cid = f"{scan}_test{n:02d}"
                    picked.append({**c, "crop_id": cid, "specimen": spec, "scan": scan, "crop": f"test{n:02d}",
                                   "size": SIZE, "role": "test_internal" if spec == "EFR" else "test_external",
                                   "test_class": cls, "touches_specimen_boundary": cls == "boundary",
                                   "z_tier": tier, "angular_half": half, "specimen_fraction_sub4": round(f, 4),
                                   "theta_deg": round(float(np.degrees(theta)), 1),
                                   "n_admissible_in_cell": len(cands),
                                   "gt_path": str(CROPS_ROOT / cid / "labels.tif"),
                                   "gray_crop_path": str(CROPS_ROOT / cid / "gray.tif")})
                    print(f"[pick] {cid} {cls:8s} tier {tier} x{c['x']} y{c['y']} z{c['z']} f={f:.3f} "
                          f"theta={np.degrees(theta):.0f} ({len(cands)} admissible)", flush=True)
        chosen += picked
        if spec == "EFR":
            dev_f = {c["crop_id"]: round(frac(geo[c["scan"]], c), 4) for c in dev["crops"]}
        geo_out.update({s: {k: v for k, v in g.items() if not k.startswith("_")} for s, g in geo.items()})
        del geo

    out = {
        "schema_version": 2,
        "description": "Phase 1f TEST set: 52 new 200^3 crops (EFR internal 6/scan; EFS, EMR, EMS, EOR, EOS external "
                       "4/scan), selected by the GT-free rule below before any cutting or evaluation. Method under "
                       "test: batos-paper-v2-fixflat (configs/batos_paper_v2_fixflat.yaml, frozen).",
        "created": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "generator": "scripts/define_test_crops.py",
        "seed": args.seed,
        "rule": __doc__.split("Rule (also stored verbatim in the JSON):")[1].strip(),
        "parameters": {"subsample": SUB, "size": SIZE, "lattice": LATTICE, "f_interior_min": F_INTERIOR,
                       "f_boundary_range": list(F_BOUNDARY), "radius_tolerance": RADIUS_TOL,
                       "z_tiers_per_scan": SPECIMENS, "per_specimen": params_spec},
        "gt_used": False,
        "geometry": geo_out,
        "dev_crops_specimen_fraction_sub4": dev_f,
        "source_volumes": sources,
        "crop_convention": dev["crop_convention"],
        "crops": chosen,
    }
    args.out.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"[done] {args.out}: {len(chosen)} crops")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
