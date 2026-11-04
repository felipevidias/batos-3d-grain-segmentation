#!/usr/bin/env python3
"""Extract only the EFRGP01_00 / EFRGP01_01 full-volume GT labels from EFR_labels.7z.

Phase 1b data preparation. Reads the archive with py7zr and extracts only the requested
members (the archive holds 55 volumes, ~2.4 GB compressed). Then reports shape, dtype and
size for each extracted volume from the TIFF header, without loading the volume.

Fails loudly if the archive or a requested member is missing. Refuses to overwrite an
existing output unless --force is given.

Example (Windows):
  C:\\BA_TOS_DATA\\envs\\batos-py312\\Scripts\\python.exe -I scripts\\extract_gt_labels.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path, PurePosixPath

DEFAULT_ARCHIVE = Path(r"C:\Users\USER\Downloads\PhDUpload\PhDUpload\EFR\Label\EFR_labels.7z")
DEFAULT_OUT = Path(r"C:\BA_TOS_DATA\gt_full")
DEFAULT_MEMBERS = ("EFRGP01_00.tif", "EFRGP01_01.tif")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--members", nargs="+", default=list(DEFAULT_MEMBERS),
                        help="Basenames of archive members to extract.")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing extracted files.")
    parser.add_argument("--sha256", action="store_true",
                        help="Also compute the SHA-256 of each extracted file (slow).")
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tiff_info(path: Path) -> dict[str, object]:
    import tifffile

    with tifffile.TiffFile(path) as handle:
        series = handle.series[0]
        return {
            "shape": [int(v) for v in series.shape],
            "dtype": str(series.dtype),
            "axes": series.axes,
            "pages": len(handle.pages),
            "bigtiff": bool(handle.is_bigtiff),
        }


def main() -> int:
    import py7zr

    args = parse_args()
    archive = args.archive.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    if not archive.is_file():
        raise FileNotFoundError(f"Archive not found: {archive}")
    out_dir.mkdir(parents=True, exist_ok=True)

    with py7zr.SevenZipFile(archive, mode="r") as handle:
        names = [n for n in handle.getnames()]
    by_base: dict[str, list[str]] = {}
    for name in names:
        by_base.setdefault(PurePosixPath(name.replace("\\", "/")).name, []).append(name)

    targets: dict[str, str] = {}
    for member in args.members:
        hits = by_base.get(member, [])
        if len(hits) != 1:
            raise RuntimeError(
                f"Expected exactly one archive member named {member}; found {hits}. "
                f"Archive has {len(names)} entries."
            )
        targets[member] = hits[0]

    for member, archive_name in targets.items():
        final = out_dir / member
        if final.exists() and not args.force:
            print(f"[skip] already extracted: {final} (use --force to overwrite)")
            continue
        start = time.time()
        print(f"[extract] {archive_name} -> {out_dir}")
        with py7zr.SevenZipFile(archive, mode="r") as handle:
            handle.extract(path=out_dir, targets=[archive_name])
        extracted = out_dir / Path(archive_name)
        if extracted != final:
            if not extracted.is_file():
                raise FileNotFoundError(f"Extraction produced no file at {extracted}")
            if final.exists():
                final.unlink()
            extracted.replace(final)
        print(f"[extract] done in {time.time() - start:.0f}s")

    report = {"archive": str(archive), "out_dir": str(out_dir), "volumes": {}}
    for member in targets:
        final = out_dir / member
        if not final.is_file():
            raise FileNotFoundError(final)
        info = tiff_info(final)
        info["bytes"] = final.stat().st_size
        if args.sha256:
            info["file_sha256"] = file_sha256(final)
        report["volumes"][member] = info
        print(f"[gt] {final}: shape={info['shape']} dtype={info['dtype']} "
              f"axes={info['axes']} bytes={info['bytes']}")

    report_path = out_dir / "extract_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"[gt] report: {report_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[extract_gt_labels][ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
