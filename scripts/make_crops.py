#!/usr/bin/env python3
"""Cut the 12 SIBGRAPI greyscale and GT crops with memory-mapped reads.

Reads data/metadata/crops_sibgrapi.json and never loads a full volume (tifffile.memmap).
Writes C:\\BA_TOS_DATA\\crops\\<crop_id>\\{gray.tif,labels.tif} and
data/metadata/crops_manifest.csv (SHA-256, shape, dtype, GT label counts, valid GT count).

The axis order comes from crop_convention.axis_order in the JSON (e.g. "zyx", meaning
volume axes 0,1,2 are indexed by the z, y, x origins). The script refuses to cut unless
crop_convention.slicing_status == "verified" (CLAUDE.md rule 3: fail loudly).

--probe-axis-orders writes nothing. It cuts the GT under each of the six orderings and
prints the label and valid-GT counts per crop, so they can be compared with historical
counts to settle the axis order.

Valid GT = label with >= 1000 voxels that does not touch any of the six crop faces.

Example (Windows):
  C:\\BA_TOS_DATA\\envs\\batos-py312\\Scripts\\python.exe -I scripts\\make_crops.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import tifffile

REPO = Path(__file__).resolve().parents[1]
DEFAULT_META = REPO / "data" / "metadata" / "crops_sibgrapi.json"
DEFAULT_MANIFEST = REPO / "data" / "metadata" / "crops_manifest.csv"
DEFAULT_OUT = Path(r"C:\BA_TOS_DATA\crops")
DEFAULT_GT_DIR = Path(r"C:\BA_TOS_DATA\gt_full")
MIN_VALID_GT = 1000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--metadata", type=Path, default=DEFAULT_META)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--gt-dir", type=Path, default=DEFAULT_GT_DIR,
                        help="Folder with full-volume GT EFRGP01_00.tif / EFRGP01_01.tif.")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--crop", action="append", default=[],
                        help="Restrict to these crop_ids (repeatable).")
    parser.add_argument("--force", action="store_true", help="Overwrite existing crops.")
    parser.add_argument("--probe-axis-orders", action="store_true",
                        help="Write nothing; report GT counts under all six orderings.")
    return parser.parse_args()


def array_sha256(array: np.ndarray) -> str:
    """SHA-256 of the C-order little-endian pixel bytes (independent of TIFF encoding)."""
    data = np.ascontiguousarray(array)
    data = data.astype(data.dtype.newbyteorder("<"), copy=False)
    return hashlib.sha256(data.tobytes()).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_memmap(path: Path) -> np.memmap:
    if not path.is_file():
        raise FileNotFoundError(path)
    volume = tifffile.memmap(str(path), mode="r")
    if volume.ndim != 3:
        raise RuntimeError(f"Expected a 3D volume, got shape {volume.shape}: {path}")
    return volume


def window(volume: np.ndarray, origin: dict[str, int], size: int, axis_order: str) -> np.ndarray:
    if sorted(axis_order) != ["x", "y", "z"]:
        raise ValueError(f"Invalid axis_order {axis_order!r}")
    starts = [int(origin[axis]) for axis in axis_order]
    for dim, start in enumerate(starts):
        if start < 0 or start + size > volume.shape[dim]:
            raise ValueError(
                f"Window {axis_order}={starts} size {size} exceeds volume shape {volume.shape}"
            )
    return np.array(volume[starts[0]:starts[0] + size,
                           starts[1]:starts[1] + size,
                           starts[2]:starts[2] + size])


def gt_stats(labels: np.ndarray) -> dict[str, int]:
    ids, counts = np.unique(labels, return_counts=True)
    keep = ids != 0
    ids, counts = ids[keep], counts[keep]
    faces = np.concatenate([
        labels[0].ravel(), labels[-1].ravel(),
        labels[:, 0].ravel(), labels[:, -1].ravel(),
        labels[:, :, 0].ravel(), labels[:, :, -1].ravel(),
    ])
    border = set(np.unique(faces[faces != 0]).tolist())
    large = counts >= MIN_VALID_GT
    touching = np.isin(ids, list(border)) if border else np.zeros(len(ids), dtype=bool)
    return {
        "n_gt_labels": int(len(ids)),
        "n_gt_touching_face": int(touching.sum()),
        "n_gt_small_lt1000": int((~large).sum()),
        "n_gt_valid": int((large & ~touching).sum()),
    }


def load_metadata(path: Path) -> dict:
    meta = json.loads(path.read_text(encoding="utf-8"))
    if len(meta.get("crops", [])) != 12:
        raise RuntimeError(f"Expected 12 crops in {path}, got {len(meta.get('crops', []))}")
    return meta


def gt_volume_path(gt_dir: Path, scan: str) -> Path:
    return gt_dir / f"{scan}.tif"


def probe(meta: dict, crops: list[dict], gt_dir: Path) -> int:
    cache: dict[str, np.ndarray] = {}
    print("crop_id," + ",".join(f"{''.join(p)}_labels/valid" for p in itertools.permutations("zyx")))
    for crop in crops:
        scan = crop["scan"]
        if scan not in cache:
            cache[scan] = open_memmap(gt_volume_path(gt_dir, scan))
        cells = []
        for perm in itertools.permutations("zyx"):
            order = "".join(perm)
            try:
                stats = gt_stats(window(cache[scan], crop, int(crop["size"]), order))
                cells.append(f"{stats['n_gt_labels']}/{stats['n_gt_valid']}")
            except ValueError:
                cells.append("out_of_bounds")
        print(crop["crop_id"] + "," + ",".join(cells), flush=True)
    return 0


def main() -> int:
    args = parse_args()
    meta = load_metadata(args.metadata.expanduser().resolve())
    crops = meta["crops"]
    if args.crop:
        wanted = set(args.crop)
        crops = [c for c in crops if c["crop_id"] in wanted]
        missing = wanted - {c["crop_id"] for c in crops}
        if missing:
            raise KeyError(f"Unknown crop_id(s): {sorted(missing)}")

    gt_dir = args.gt_dir.expanduser().resolve()
    if args.probe_axis_orders:
        return probe(meta, crops, gt_dir)

    convention = meta.get("crop_convention", {})
    if convention.get("slicing_status") != "verified":
        raise RuntimeError(
            "crop_convention.slicing_status is not 'verified' in crops_sibgrapi.json. "
            "Settle the axis order (Q7) before cutting crops."
        )
    axis_order = convention["axis_order"]

    out_root = args.out_root.expanduser().resolve()
    gray_cache: dict[str, np.ndarray] = {}
    gt_cache: dict[str, np.ndarray] = {}
    rows = []
    for crop in crops:
        crop_id, scan, size = crop["crop_id"], crop["scan"], int(crop["size"])
        out_dir = out_root / crop_id
        gray_path, label_path = out_dir / "gray.tif", out_dir / "labels.tif"
        if (gray_path.exists() or label_path.exists()) and not args.force:
            raise FileExistsError(f"{out_dir} already has crops; use --force to overwrite")
        if scan not in gray_cache:
            gray_cache[scan] = open_memmap(Path(meta["source_volumes"][scan]["gray_path"]))
            gt_cache[scan] = open_memmap(gt_volume_path(gt_dir, scan))
        gray = window(gray_cache[scan], crop, size, axis_order)
        labels = window(gt_cache[scan], crop, size, axis_order)
        out_dir.mkdir(parents=True, exist_ok=True)
        tifffile.imwrite(gray_path, gray)
        tifffile.imwrite(label_path, labels)
        stats = gt_stats(labels)
        row = {
            "crop_id": crop_id, "scan": scan, "role": crop.get("role", ""),
            "x": crop["x"], "y": crop["y"], "z": crop["z"], "size": size,
            "axis_order": axis_order,
            "gray_path": str(gray_path), "gray_shape": "x".join(map(str, gray.shape)),
            "gray_dtype": str(gray.dtype), "gray_pixel_sha256": array_sha256(gray),
            "gray_file_sha256": file_sha256(gray_path),
            "gray_min": int(gray.min()), "gray_max": int(gray.max()),
            "labels_path": str(label_path), "labels_shape": "x".join(map(str, labels.shape)),
            "labels_dtype": str(labels.dtype), "labels_pixel_sha256": array_sha256(labels),
            "labels_file_sha256": file_sha256(label_path),
            **stats,
        }
        rows.append(row)
        print(f"[crop] {crop_id}: gray {row['gray_shape']} {row['gray_dtype']}, "
              f"labels {row['labels_dtype']}, GT labels={stats['n_gt_labels']}, "
              f"valid={stats['n_gt_valid']}", flush=True)

    manifest = args.manifest.expanduser().resolve()
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    total_valid = sum(r["n_gt_valid"] for r in rows)
    print(f"[crop] manifest: {manifest}")
    print(f"[crop] total valid GT over {len(rows)} crops: {total_valid}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[make_crops][ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
