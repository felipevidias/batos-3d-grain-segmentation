#!/usr/bin/env python3
"""Phase 5i (only after the v3c freeze): cut TEST SET 2 (data/metadata/crops_test2.json) for v3 and v3c.

Per specimen:
  1. extract its two label volumes once (EFR: the kept C:\\BA_TOS_DATA\\gt_full copies);
  2. cut the 200^3 probes with the unchanged test-2 cutter (scripts/cut_test2.py -> extract_and_cut.py), which also
     writes the standard manifest data/metadata/crops_manifest_test2.csv;            -> v3 runs on these
  3. per crop: r_probe = batos_v3.scale.estimate(probe grey) (GT-free); native side S = round(200 r_probe / 10);
     window origin = probe origin + floor((200 - S) / 2) per axis (same centre; exact when S is even), shifted by the
     minimum amount needed to lie inside the volume (shift recorded; S larger than the volume -> fail loudly);
     grey + GT cut from the full volumes                                                -> v3c runs on these
     valid GT for v3c: >= round(1000 (r_probe / 10)^3) native voxels, not touching a face;
  4. delete the extracted label volumes of the specimen (disk budget).
Outputs: C:\\BA_TOS_DATA\\crops_test2\\<id> (probes), C:\\BA_TOS_DATA\\crops_test2_v3c\\<id> (windows),
data/metadata/crops_test2_v3c.json, data/metadata/crops_manifest_test2_v3c.csv. Restartable (done crops skipped).
"""

from __future__ import annotations

import csv
import gc
import hashlib
import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import tifffile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v3.scale import estimate  # noqa: E402
from batos_v3c import eval_lean  # noqa: E402

_spec = importlib.util.spec_from_file_location("cut_test2", REPO / "scripts" / "cut_test2.py")
ct2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ct2)
eac = ct2.eac  # extract_and_cut re-pointed at test set 2

META = REPO / "data" / "metadata" / "crops_test2.json"
OUT_JSON = REPO / "data" / "metadata" / "crops_test2_v3c.json"
OUT_CSV = REPO / "data" / "metadata" / "crops_manifest_test2_v3c.csv"
ROOT = Path(r"C:\BA_TOS_DATA\crops_test2_v3c")
PROBE, R_REF = 200, 10.0


def sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def main() -> int:
    meta = json.loads(META.read_text(encoding="utf-8"))
    order = meta["crop_convention"]["axis_order"]
    if order != "zyx":
        raise RuntimeError(f"axis order {order}")
    done = json.loads(OUT_JSON.read_text(encoding="utf-8"))["crops"] if OUT_JSON.exists() else []
    have = {c["crop_id"] for c in done}
    specimens = list(dict.fromkeys(c["specimen"] for c in meta["crops"]))
    for spec in specimens:
        crops = [c for c in meta["crops"] if c["specimen"] == spec]
        todo = [c for c in crops if c["crop_id"] not in have]
        if not todo:
            print(f"[skip] {spec}: all v3c windows cut", flush=True)
            continue
        scans = sorted({c["scan"] for c in crops})
        if spec == "EFR":
            gt_paths, extracted = {s: eac.EFR_GT / f"{s}.tif" for s in scans}, []
        else:
            ext = eac.extract(Path(meta["source_volumes"][scans[0]]["labels_archive"]), [f"{s}.tif" for s in scans],
                              eac.TMP / spec)
            gt_paths, extracted = {s: ext[f"{s}.tif"] for s in scans}, list(ext.values())
        try:
            argv = sys.argv
            sys.argv = ["cut_test2.py", "--specimens", spec, "--keep-extracted"]
            try:
                eac.main()  # probes + standard manifest (skips crops already in it)
            finally:
                sys.argv = argv
            for c in todo:
                probe = tifffile.imread(c["gray_crop_path"])
                if probe.shape != (PROBE,) * 3:
                    raise RuntimeError(f"{c['crop_id']}: probe shape {probe.shape}")
                r = float(estimate(probe)[0]["r_granulo"])
                del probe
                S = int(round(PROBE * r / R_REF))
                gray_v = eac.mc.open_memmap(Path(meta["source_volumes"][c["scan"]]["gray_path"]))
                gt_v = eac.mc.open_memmap(gt_paths[c["scan"]])
                if gray_v.shape != gt_v.shape:
                    raise RuntimeError(f"{c['crop_id']}: gray {gray_v.shape} != labels {gt_v.shape}")
                origin, shift = {}, {}
                for ax, n in zip("zyx", gray_v.shape):
                    if S > n:
                        raise RuntimeError(f"{c['crop_id']}: native side {S} > volume extent {n} on {ax}")
                    o = int(c[ax]) + (PROBE - S) // 2
                    origin[ax] = min(max(o, 0), n - S)
                    shift[ax] = origin[ax] - o
                gray = eac.mc.window(gray_v, origin, S, order)
                lab = eac.mc.window(gt_v, origin, S, order)
                d = ROOT / c["crop_id"]
                d.mkdir(parents=True, exist_ok=True)
                tifffile.imwrite(d / "gray.tif", gray)
                tifffile.imwrite(d / "labels.tif", lab)
                mgv = int(round(1000 * (r / R_REF) ** 3))
                n_valid = int(len(eval_lean.protocol(lab, mgv)["valid_ids"]))
                rec = {"crop_id": c["crop_id"], "specimen": spec, "scan": c["scan"], "test_class": c["test_class"],
                       "probe_origin_zyx": [int(c[a]) for a in "zyx"], "r_probe": r, "native_side": S,
                       "f_canon": R_REF / r, "origin_zyx": [origin[a] for a in "zyx"],
                       "shift_zyx": [shift[a] for a in "zyx"], "centred_exactly": (PROBE - S) % 2 == 0 and
                       not any(shift.values()), "min_gt_voxels": mgv, "valid_gt": n_valid,
                       "gray_crop_path": str(d / "gray.tif"), "gt_path": str(d / "labels.tif"),
                       "probe_gray_path": c["gray_crop_path"], "probe_gt_path": c["gt_path"],
                       "gray_sha256": sha(gray), "labels_sha256": sha(lab)}
                done.append(rec)
                have.add(c["crop_id"])
                print(f"[v3c] {c['crop_id']}: r_probe {r:g} -> native {S}^3 origin {rec['origin_zyx']} shift "
                      f"{rec['shift_zyx']}, valid GT {n_valid} (min {mgv})", flush=True)
                del gray, lab, gray_v, gt_v
                OUT_JSON.write_text(json.dumps({"schema_version": 1, "source": str(META.relative_to(REPO)),
                                                "description": __doc__, "crops": done}, indent=1), encoding="utf-8")
        finally:
            if extracted:
                gc.collect()
                for p in extracted:
                    p.unlink(missing_ok=True)
                shutil.rmtree(eac.TMP / spec, ignore_errors=True)
                print(f"[clean] deleted {[p.name for p in extracted]}", flush=True)
    cols = ["crop_id", "specimen", "scan", "test_class", "r_probe", "native_side", "origin_zyx", "shift_zyx",
            "centred_exactly", "min_gt_voxels", "valid_gt", "gray_sha256", "labels_sha256"]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(done)
    print(f"[done] {len(done)} v3c windows -> {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
