#!/usr/bin/env python3
"""Define held-out TEST SET 2 by the same GT-free rule as test set 1 and write data/metadata/crops_test2.json
(seed 20261109). Written in Phase 5f step 4 (never run there); run in Phase 5i after the batos-paper-v3c freeze.

The rule is that of scripts/define_test_crops.py (its geometry / fraction / overlap functions and constants are imported
unchanged; rule text stored verbatim in the JSON), with ONE addition to step 4: a box may not overlap any test-set-1 box
of the same specimen (data/metadata/crops_test.json), in addition to already-selected set-2 boxes and, for EFR, the 6
development origins. Same strata (z-tier x angular half x class), same cell order, one numpy default_rng(SEED).
Crop ids: <scan>_t2_NN. Crops: C:\\BA_TOS_DATA\\crops_test2\\<crop_id>\\{gray,labels}.tif (cut by scripts/cut_test2.py).
Reads only greyscale volumes (subsampled); never GT.
"""

from __future__ import annotations

import argparse
import datetime as dt
import importlib.util
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("dtc", REPO / "scripts" / "define_test_crops.py")
dtc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dtc)

TEST1 = REPO / "data" / "metadata" / "crops_test.json"
OUT = REPO / "data" / "metadata" / "crops_test2.json"
CROPS_ROOT = Path(r"C:\BA_TOS_DATA\crops_test2")
SEED = 20261109


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    if args.out.exists():
        raise FileExistsError(f"{args.out} exists; the test set is defined once")
    dev = json.loads(dtc.DEV.read_text(encoding="utf-8"))
    t1 = json.loads(TEST1.read_text(encoding="utf-8"))
    dev_boxes = [dict(x=x, y=y, z=z) for x, y, z in sorted({(c["x"], c["y"], c["z"]) for c in dev["crops"]})]
    rng = np.random.default_rng(args.seed)
    chosen, geo_out, params_spec = [], {}, {}
    for spec, n_tiers in dtc.SPECIMENS.items():
        scans = {0: f"{spec}GP01_00", 1: f"{spec}GP01_01"}
        geo = {}
        for scan in scans.values():
            geo[scan] = dtc.geometry(t1["source_volumes"][scan]["gray_path"])
            g = geo[scan]
            if g["mask_sha256"] != t1["geometry"][scan]["mask_sha256"]:
                raise RuntimeError(f"{scan}: specimen mask differs from test set 1")
            print(f"[geo] {scan}: z_ok {g['z_ok']}, axis {g['axis_xy'][0]:.1f},{g['axis_xy'][1]:.1f} "
                  f"(mask identical to test set 1)", flush=True)
        shape = geo[scans[0]]["shape"]
        z_lo = max(g["z_ok"][0] for g in geo.values())
        z_hi = min(g["z_ok"][1] for g in geo.values())
        z_min = int(np.ceil(z_lo / dtc.LATTICE) * dtc.LATTICE)
        z_max = int(np.floor((z_hi + dtc.SUB - dtc.SIZE) / dtc.LATTICE) * dtc.LATTICE)
        tiers = np.linspace(z_min, z_max + dtc.LATTICE, n_tiers + 1)
        t1_boxes = [dict(x=c["x"], y=c["y"], z=c["z"]) for c in t1["crops"] if c["specimen"] == spec]
        excluded = (dev_boxes if spec == "EFR" else []) + t1_boxes
        params_spec[spec] = {"admissible_origin_z": [z_min, z_max], "z_tier_edges": tiers.tolist(),
                             "scan_of_angular_half": scans,
                             "excluded_boxes": ("12 development crops (6 origins) + " if spec == "EFR" else "")
                             + f"{len(t1_boxes)} test-set-1 boxes"}
        lattice = [dict(x=x, y=y, z=z) for z in range(z_min, z_max + 1, dtc.LATTICE)
                   for y in range(0, shape[1] - dtc.SIZE + 1, dtc.LATTICE)
                   for x in range(0, shape[2] - dtc.SIZE + 1, dtc.LATTICE)]
        picked = []
        for half, scan in scans.items():
            ax, ay = geo[scan]["axis_xy"]
            for tier in range(n_tiers):
                for cls in dtc.CLASSES:
                    cands = []
                    for c in lattice:
                        if not (tiers[tier] <= c["z"] < tiers[tier + 1]):
                            continue
                        theta = np.arctan2(c["y"] + dtc.SIZE / 2 - ay, c["x"] + dtc.SIZE / 2 - ax)
                        if (0 if theta >= 0 else 1) != half:
                            continue
                        if any(dtc.overlaps(c, b) for b in excluded + picked):
                            continue
                        f = dtc.frac(geo[scan], c)
                        if (cls == "interior" and f >= dtc.F_INTERIOR) or \
                                (cls == "boundary" and dtc.F_BOUNDARY[0] <= f <= dtc.F_BOUNDARY[1]):
                            cands.append((c, f, theta))
                    if not cands:
                        raise RuntimeError(f"no admissible candidate for {scan} tier {tier} {cls}")
                    c, f, theta = cands[int(rng.integers(len(cands)))]
                    n = len(picked) + 1
                    cid = f"{scan}_t2_{n:02d}"
                    picked.append({**c, "crop_id": cid, "specimen": spec, "scan": scan, "crop": f"t2_{n:02d}",
                                   "size": dtc.SIZE, "role": "test2_internal" if spec == "EFR" else "test2_external",
                                   "test_class": cls, "touches_specimen_boundary": cls == "boundary",
                                   "z_tier": tier, "angular_half": half, "specimen_fraction_sub4": round(f, 4),
                                   "theta_deg": round(float(np.degrees(theta)), 1), "n_admissible_in_cell": len(cands),
                                   "gt_path": str(CROPS_ROOT / cid / "labels.tif"),
                                   "gray_crop_path": str(CROPS_ROOT / cid / "gray.tif")})
                    print(f"[pick] {cid} {cls:8s} tier {tier} x{c['x']} y{c['y']} z{c['z']} f={f:.3f} "
                          f"({len(cands)} admissible)", flush=True)
        chosen += picked
        geo_out.update({s: {k: v for k, v in g.items() if not k.startswith("_")} for s, g in geo.items()})
        del geo
    for a in chosen:  # final guard: no overlap with dev (EFR) or any test-set-1 box of the same specimen
        bad = [b["crop_id"] for b in t1["crops"] if b["specimen"] == a["specimen"] and dtc.overlaps(a, b)]
        if bad:
            raise RuntimeError(f"{a['crop_id']} overlaps test set 1: {bad}")
    out = {
        "schema_version": 2,
        "description": "TEST SET 2: 52 new 200^3 crops (EFR 6/scan; EFS, EMR, EMS, EOR, EOS 4/scan), same GT-free "
                       "rule as test set 1 plus no overlap with any test-set-1 box. Defined in Phase 5i after the "
                       "batos-paper-v3c freeze (v3.1 of Phase 5f was never frozen). Methods: v3 (200^3 probe) and v3c "
                       "(re-cut window of native side round(200 r/10), scripts/cut_test2_v3c.py).",
        "created": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "generator": "scripts/define_test_crops2.py", "seed": args.seed,
        "rule": t1["rule"] + "\n 8. (test set 2) A box may also not overlap any test-set-1 box of the same specimen "
                             "(data/metadata/crops_test.json).",
        "parameters": {**t1["parameters"], "per_specimen": params_spec}, "gt_used": False, "geometry": geo_out,
        "source_volumes": t1["source_volumes"], "crop_convention": t1["crop_convention"], "crops": chosen,
    }
    args.out.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"[done] {args.out}: {len(chosen)} crops")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
