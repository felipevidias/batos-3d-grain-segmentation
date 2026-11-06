#!/usr/bin/env python3
"""Phase 1f: extract the _00/_01 label volumes, cut the 52 test crops, write the test manifest, delete the extracts.

For each specimen in data/metadata/crops_test.json:
  1. extract ONLY <S>GP01_00.tif and <S>GP01_01.tif from <S>_labels.7z with py7zr into C:\\BA_TOS_DATA\\gt_full_tmp\\
     (EFR is not extracted: the Phase 1b copies in C:\\BA_TOS_DATA\\gt_full\\ are used and kept as they are);
  2. cut gray + GT 200^3 crops by memmap (vol[z:z+200, y:y+200, x:x+200], zyx as crops_sibgrapi.json) into
     C:\\BA_TOS_DATA\\crops_test\\<crop_id>\\{gray,labels}.tif;
  3. delete the extracted full label volumes of that specimen before the next one (disk budget).
Manifest data/metadata/crops_manifest_test.csv: pixel/file SHA-256, shape, dtype, GT label counts, valid GT
(>= 1000 voxels, not touching a face), e13 Otsu foreground fraction (GT-free), GT-labelled fraction, and whether the
crop intersects the specimen boundary (by the selection rule, not by GT).
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import py7zr
import tifffile

REPO = Path(__file__).resolve().parents[1]
META = REPO / "data" / "metadata" / "crops_test.json"
MANIFEST = REPO / "data" / "metadata" / "crops_manifest_test.csv"
TMP = Path(r"C:\BA_TOS_DATA\gt_full_tmp")
EFR_GT = Path(r"C:\BA_TOS_DATA\gt_full")
OUT = Path(r"C:\BA_TOS_DATA\crops_test")

spec = importlib.util.spec_from_file_location("make_crops", REPO / "scripts" / "make_crops.py")
mc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mc)
sys.path.insert(0, str(REPO / "src"))
from batos_v2.scale_granulo import load_e13  # noqa: E402  (e13 foreground, imported read-only)


def extract(archive: Path, members: list[str], out_dir: Path) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {m: out_dir / m for m in members}
    todo = [m for m, p in paths.items() if not p.exists()]
    if todo:
        t0 = time.time()
        with py7zr.SevenZipFile(archive, mode="r") as z:
            names = set(z.getnames())
            missing = set(todo) - names
            if missing:
                raise FileNotFoundError(f"{archive}: members not found {sorted(missing)}")
            z.extract(path=out_dir, targets=todo)
        print(f"[extract] {archive.name}: {todo} in {time.time() - t0:.0f}s", flush=True)
    for m, p in paths.items():
        if not p.is_file():
            raise FileNotFoundError(p)
    return paths


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--specimens", nargs="*", default=None)
    ap.add_argument("--keep-extracted", action="store_true", help="do not delete gt_full_tmp volumes (debug only)")
    args = ap.parse_args()
    meta = json.loads(META.read_text(encoding="utf-8"))
    order = meta["crop_convention"]["axis_order"]
    if meta["crop_convention"].get("slicing_status") != "verified":
        raise RuntimeError("axis order not verified")
    e13 = load_e13()
    rows = []
    if MANIFEST.exists():
        rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8")))
    done = {r["crop_id"] for r in rows}
    specimens = list(dict.fromkeys(c["specimen"] for c in meta["crops"]))
    for specimen in args.specimens or specimens:
        crops = [c for c in meta["crops"] if c["specimen"] == specimen and c["crop_id"] not in done]
        if not crops:
            print(f"[skip] {specimen}: all crops already in the manifest", flush=True)
            continue
        scans = sorted({c["scan"] for c in crops})
        if specimen == "EFR":
            gt_paths = {s: EFR_GT / f"{s}.tif" for s in scans}
            extracted = []
        else:
            archive = Path(meta["source_volumes"][scans[0]]["labels_archive"])
            ext = extract(archive, [f"{s}.tif" for s in scans], TMP / specimen)
            gt_paths = {s: ext[f"{s}.tif"] for s in scans}
            extracted = list(ext.values())
        try:
            for c in crops:
                gray_v = mc.open_memmap(Path(meta["source_volumes"][c["scan"]]["gray_path"]))
                gt_v = mc.open_memmap(gt_paths[c["scan"]])
                if gray_v.shape != gt_v.shape:
                    raise RuntimeError(f"{c['crop_id']}: gray {gray_v.shape} != labels {gt_v.shape}")
                gray = mc.window(gray_v, c, c["size"], order)
                labels = mc.window(gt_v, c, c["size"], order)
                d = OUT / c["crop_id"]
                d.mkdir(parents=True, exist_ok=True)
                gp, lp = d / "gray.tif", d / "labels.tif"
                tifffile.imwrite(gp, gray)
                tifffile.imwrite(lp, labels)
                fg, _ = e13.make_foreground(e13.robust_u8(gray))
                stats = mc.gt_stats(labels)
                row = {"crop_id": c["crop_id"], "specimen": specimen, "scan": c["scan"], "role": c["role"],
                       "test_class": c["test_class"], "x": c["x"], "y": c["y"], "z": c["z"], "size": c["size"],
                       "axis_order": order, "gray_path": str(gp), "gray_shape": "x".join(map(str, gray.shape)),
                       "gray_dtype": str(gray.dtype), "gray_pixel_sha256": mc.array_sha256(gray),
                       "gray_file_sha256": mc.file_sha256(gp), "gray_min": int(gray.min()),
                       "gray_max": int(gray.max()), "labels_path": str(lp),
                       "labels_shape": "x".join(map(str, labels.shape)), "labels_dtype": str(labels.dtype),
                       "labels_pixel_sha256": mc.array_sha256(labels), "labels_file_sha256": mc.file_sha256(lp),
                       "labels_source": str(gt_paths[c["scan"]]), **stats,
                       "foreground_fraction_e13_otsu": round(float(fg.mean()), 4),
                       "gt_labelled_fraction": round(float((labels > 0).mean()), 4),
                       "specimen_fraction_sub4": c["specimen_fraction_sub4"],
                       "touches_specimen_boundary": c["touches_specimen_boundary"]}
                rows.append(row)
                print(f"[crop] {c['crop_id']} {c['test_class']}: labels {stats['n_gt_labels']}, "
                      f"valid {stats['n_gt_valid']}, fg {row['foreground_fraction_e13_otsu']:.3f}, "
                      f"gt {row['gt_labelled_fraction']:.3f}", flush=True)
                with MANIFEST.open("w", encoding="utf-8", newline="") as h:
                    w = csv.DictWriter(h, fieldnames=list(rows[-1].keys()))
                    w.writeheader()
                    w.writerows(rows)
                del gray_v, gt_v
        finally:
            if extracted and not args.keep_extracted:
                import gc
                gc.collect()
                for p in extracted:
                    p.unlink(missing_ok=True)
                shutil.rmtree(TMP / specimen, ignore_errors=True)
                print(f"[clean] deleted {[p.name for p in extracted]}", flush=True)
    print(f"[done] {MANIFEST}: {len(rows)} crops, valid GT total {sum(int(r['n_gt_valid']) for r in rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
