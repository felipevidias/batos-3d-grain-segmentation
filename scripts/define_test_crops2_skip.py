#!/usr/bin/env python3
"""Define TEST SET 2 (seed 20261109) with the rule of scripts/define_test_crops2.py plus ONE fallback, fixed before any
outcome was seen: a stratum cell (scan x z-tier x class) with no admissible candidate is SKIPPED and recorded (no
relaxation of any constraint, no random draw for that cell); everything else is identical.

Why: the frozen define_test_crops2.py (hashed in configs/batos_paper_v3c.yaml, therefore not edited) raised
"no admissible candidate for EFRGP01_00 tier 1 interior" on 2026-11-10 -- after the dev 6 origins, the 6+6 test-set-1
EFR boxes and the newly picked set-2 boxes are excluded, the small EFR interior region is exhausted. It wrote nothing.
This script reuses the frozen module's constants and the define_test_crops.py geometry / fraction / overlap functions
unchanged, and writes the same file (data/metadata/crops_test2.json), which the frozen cut / run / aggregate drivers read.
Reads only greyscale volumes (subsampled); never GT.
"""

from __future__ import annotations

import datetime as dt
import importlib.util
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("dtc2", REPO / "scripts" / "define_test_crops2.py")
dtc2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dtc2)
dtc = dtc2.dtc

RULE_FALLBACK = ("\n 9. (test set 2, fallback fixed before any outcome, 2026-11-10) A stratum cell with no admissible "
                 "candidate is skipped and recorded in 'skipped_cells'; no constraint is relaxed and no random number is "
                 "drawn for it.")


def main() -> int:
    out_path = dtc2.OUT
    if out_path.exists():
        raise FileExistsError(f"{out_path} exists; the test set is defined once")
    dev = json.loads(dtc.DEV.read_text(encoding="utf-8"))
    t1 = json.loads(dtc2.TEST1.read_text(encoding="utf-8"))
    dev_boxes = [dict(x=x, y=y, z=z) for x, y, z in sorted({(c["x"], c["y"], c["z"]) for c in dev["crops"]})]
    rng = np.random.default_rng(dtc2.SEED)
    chosen, skipped, geo_out, params_spec = [], [], {}, {}
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
                        skipped.append({"specimen": spec, "scan": scan, "z_tier": tier, "angular_half": half,
                                        "test_class": cls, "reason": "no admissible candidate"})
                        print(f"[skip] {scan} tier {tier} {cls}: no admissible candidate (cell skipped)", flush=True)
                        continue
                    c, f, theta = cands[int(rng.integers(len(cands)))]
                    n = len(picked) + 1
                    cid = f"{scan}_t2_{n:02d}"
                    picked.append({**c, "crop_id": cid, "specimen": spec, "scan": scan, "crop": f"t2_{n:02d}",
                                   "size": dtc.SIZE, "role": "test2_internal" if spec == "EFR" else "test2_external",
                                   "test_class": cls, "touches_specimen_boundary": cls == "boundary",
                                   "z_tier": tier, "angular_half": half, "specimen_fraction_sub4": round(f, 4),
                                   "theta_deg": round(float(np.degrees(theta)), 1), "n_admissible_in_cell": len(cands),
                                   "gt_path": str(dtc2.CROPS_ROOT / cid / "labels.tif"),
                                   "gray_crop_path": str(dtc2.CROPS_ROOT / cid / "gray.tif")})
                    print(f"[pick] {cid} {cls:8s} tier {tier} x{c['x']} y{c['y']} z{c['z']} f={f:.3f} "
                          f"({len(cands)} admissible)", flush=True)
        chosen += picked
        geo_out.update({s: {k: v for k, v in g.items() if not k.startswith("_")} for s, g in geo.items()})
        del geo
    for a in chosen:  # final guard: no overlap with dev (EFR) or any test-set-1 box of the same specimen
        bad = [b["crop_id"] for b in t1["crops"] if b["specimen"] == a["specimen"] and dtc.overlaps(a, b)]
        if bad:
            raise RuntimeError(f"{a['crop_id']} overlaps test set 1: {bad}")
    n_spec = {s: sum(c["specimen"] == s for c in chosen) for s in dtc.SPECIMENS}
    out = {
        "schema_version": 2,
        "description": f"TEST SET 2: {len(chosen)} new 200^3 crops ({n_spec}; {len(skipped)} stratum cells skipped, "
                       "see skipped_cells), same GT-free rule as test set 1 plus no overlap with any test-set-1 box. "
                       "Defined in Phase 5i after the batos-paper-v3c freeze. Methods: v3 (200^3 probe) and v3c "
                       "(re-cut window of native side round(200 r/10), scripts/cut_test2_v3c.py).",
        "created": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "generator": "scripts/define_test_crops2_skip.py (rule of scripts/define_test_crops2.py + fallback 9)",
        "seed": dtc2.SEED,
        "rule": t1["rule"] + "\n 8. (test set 2) A box may also not overlap any test-set-1 box of the same specimen "
                             "(data/metadata/crops_test.json)." + RULE_FALLBACK,
        "parameters": {**t1["parameters"], "per_specimen": params_spec}, "gt_used": False, "geometry": geo_out,
        "source_volumes": t1["source_volumes"], "crop_convention": t1["crop_convention"], "crops": chosen,
        "skipped_cells": skipped,
    }
    out_path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"[done] {out_path}: {len(chosen)} crops, {len(skipped)} cells skipped {n_spec}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
