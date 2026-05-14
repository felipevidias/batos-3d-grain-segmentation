#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys
import csv
import pandas as pd

ROOT = Path(".")
SCRIPT_DIR = ROOT / "code/segmentation/batos"
OUTROOT = ROOT / "results_batos/paper_ablation_full"
TABLEDIR = ROOT / "results_batos/paper_tables/final_revision_clean"

OUTROOT.mkdir(parents=True, exist_ok=True)
TABLEDIR.mkdir(parents=True, exist_ok=True)

PY = sys.executable

CROPS = {
    "crop02": ("600", "800", "300"),
    "crop03": ("900", "600", "300"),
    "crop04": ("500", "500", "700"),
    "crop05": ("800", "900", "700"),
    "crop06": ("350", "650", "500"),
    "crop07": ("700", "400", "500"),
}

def run(cmd):
    print("RUN:", " ".join(map(str, cmd)))
    subprocess.run([str(x) for x in cmd], check=True)

def find_one(pattern):
    hits = sorted(ROOT.glob(pattern))
    hits = [p for p in hits if p.is_file()]
    return hits[-1] if hits else None

def crop_paths(scan, crop, x, y, z):
    if scan == "EFRGP01_00":
        if crop in ["crop02", "crop03", "crop04", "crop05"]:
            base = ROOT / f"results_batos/EFRGP01_00_multicrop/{crop}_x{x}_y{y}_z{z}_s200"
        else:
            base = ROOT / f"results_batos/EFRGP01_00_extra_crops/{crop}_x{x}_y{y}_z{z}_s200"

        gray = base / "00_input" / f"{crop}_x{x}_y{y}_z{z}_s200_gray_crop_8u.tif"
        gt = base / "00_input" / f"{crop}_x{x}_y{y}_z{z}_s200_label_crop.tif"

        full = ROOT / f"results_batos/paper_tables/batos_size_filter_sweep/min500/{crop}/{crop}_batos_min500_labels.tif"
        raw = base / "batos_pipeline" / "05_watershed" / "batos_distance_labels.tif"

    else:
        base = ROOT / f"results_batos/EFRGP01_01_selected_crops/{crop}_x{x}_y{y}_z{z}_s200"

        gray = base / "00_input" / f"EFRGP01_01_{crop}_x{x}_y{y}_z{z}_s200_gray_crop_8u.tif"
        gt = base / "00_input" / f"EFRGP01_01_{crop}_x{x}_y{y}_z{z}_s200_label_crop.tif"

        full = base / "batos_pipeline" / "05_watershed" / "batos_distance_filtered_min500_labels.tif"
        raw = base / "batos_pipeline" / "05_watershed" / "batos_distance_labels.tif"

    return base, gray, gt, full, raw

def run_tos_seed(flat, outdir, stem):
    seed_dir = outdir / "02_tos_seed"
    seed_dir.mkdir(parents=True, exist_ok=True)

    run([
        PY, SCRIPT_DIR / "batos_tos_seed.py", flat,
        "--out-dir", seed_dir,
        "--height-factor", "0.012",
        "--area-factor", "22",
        "--area-min", "2",
        "--area-max", "2200",
        "--delta-parent-min", "0",
        "--dark-percentile-max", "78",
        "--min-score", "0.35",
        "--max-selected-nodes", "900",
        "--support-area-min", "20",
        "--support-area-max", "45000",
        "--support-max-steps", "16",
        "--component-min", "2",
        "--component-max", "1500",
    ])

    raw_seed = find_one(str(seed_dir / "*_tosconop_autodual_seed_raw.tif"))
    if raw_seed is None:
        raise RuntimeError(f"raw ToS seed not found in {seed_dir}")
    return raw_seed

def run_background(original, flat, raw_seed, outdir, stem):
    bg_dir = outdir / "03_background"
    bg_dir.mkdir(parents=True, exist_ok=True)

    run([
        PY, SCRIPT_DIR / "batos_background.py", original,
        "--flattened", flat,
        "--seed", raw_seed,
        "--out-dir", bg_dir,
        "--prefix", stem,
        "--dark-percentile", "40",
        "--otsu-offset", "-10",
        "--dark-threshold-max", "135",
        "--dark-connectivity", "1",
        "--seed-area-min", "1",
        "--seed-area-max", "3000",
        "--max-external-overlap", "0.10",
        "--min-internal-overlap", "0.00",
        "--internal-search-radius", "6",
        "--center-dark-max", "210",
        "--min-ring-contrast", "-5",
        "--max-ring-external-coverage", "0.70",
        "--nms-min-dist", "5",
    ])

    internal = find_one(str(bg_dir / "*_internal_dark_candidates.tif"))
    if internal is None:
        raise RuntimeError(f"internal candidates not found in {bg_dir}")
    return internal

def run_markers(original, flat, internal, outdir, stem):
    mk_dir = outdir / "04_markers"
    mk_dir.mkdir(parents=True, exist_ok=True)

    run([
        PY, SCRIPT_DIR / "batos_markers.py", original,
        "--flattened", flat,
        "--internal-dark", internal,
        "--out-dir", mk_dir,
        "--prefix", stem,
        "--connectivity", "1",
        "--area-min", "3",
        "--area-max", "900",
        "--max-dim", "30",
        "--max-elongation", "8",
        "--min-fill", "0.005",
        "--center-dark-max", "205",
        "--min-ring-contrast", "-2",
        "--nms-min-dist", "7",
        "--max-markers", "420",
    ])

    markers = find_one(str(mk_dir / "*_batos_marker_components_raw.tif"))
    if markers is None:
        raise RuntimeError(f"marker components not found in {mk_dir}")
    return markers

def run_watershed(original, markers, outdir):
    ws_dir = outdir / "05_watershed"
    ws_dir.mkdir(parents=True, exist_ok=True)

    run([
        PY, SCRIPT_DIR / "batos_watershed_distance.py",
        "--image", original,
        "--markers", markers,
        "--out-dir", ws_dir,
        "--prefix", "batos_distance",
        "--threshold-offset", "-10",
        "--close-radius", "1",
        "--force-markers-foreground",
        "--force-radius", "1",
    ])

    raw_labels = ws_dir / "batos_distance_labels.tif"
    if not raw_labels.exists():
        raise RuntimeError(f"raw watershed labels not found: {raw_labels}")

    run([
        PY, SCRIPT_DIR / "filter_predicted_labels_by_size_labels_only.py",
        "--labels", raw_labels,
        "--out-dir", ws_dir,
        "--prefix", "batos_distance_filtered",
        "--min-size", "500",
        "--relabel",
    ])

    filtered = ws_dir / "batos_distance_filtered_min500_labels.tif"
    if not filtered.exists():
        raise RuntimeError(f"filtered labels not found: {filtered}")

    return raw_labels, filtered

def run_variant_no_boundary(scan, crop, gray):
    outdir = OUTROOT / "without_boundary_rejection" / scan / crop
    stem = f"{scan}_{crop}_no_boundary"

    if (outdir / "05_watershed/batos_distance_filtered_min500_labels.tif").exists():
        return outdir / "05_watershed/batos_distance_labels.tif", outdir / "05_watershed/batos_distance_filtered_min500_labels.tif"

    # Full flattening, full ToS, but skip batos_background.py.
    pp_dir = outdir / "01_preprocess"
    pp_dir.mkdir(parents=True, exist_ok=True)

    run([
        PY, SCRIPT_DIR / "batos_preprocess.py", gray,
        "--out-dir", pp_dir,
        "--prefix", stem,
        "--lambda-xy", "15",
    ])

    flat = find_one(str(pp_dir / "*_flatten_lambda*_8u.tif"))
    if flat is None:
        raise RuntimeError(f"flatten not found in {pp_dir}")

    raw_seed = run_tos_seed(flat, outdir, stem)

    # Ablation: treat raw ToS seed as internal-dark candidates.
    markers = run_markers(gray, flat, raw_seed, outdir, stem)
    return run_watershed(gray, markers, outdir)

def run_variant_no_flattening(scan, crop, gray):
    outdir = OUTROOT / "without_flattening" / scan / crop
    stem = f"{scan}_{crop}_no_flattening"

    if (outdir / "05_watershed/batos_distance_filtered_min500_labels.tif").exists():
        return outdir / "05_watershed/batos_distance_labels.tif", outdir / "05_watershed/batos_distance_filtered_min500_labels.tif"

    # Ablation: use original 8-bit gray crop as the flattened volume.
    flat = gray

    raw_seed = run_tos_seed(flat, outdir, stem)
    internal = run_background(gray, flat, raw_seed, outdir, stem)
    markers = run_markers(gray, flat, internal, outdir, stem)
    return run_watershed(gray, markers, outdir)

manifest = []

def add_manifest(scan, crop, variant, gt, pred):
    manifest.append({
        "scan": scan,
        "crop_id": crop,
        "variant": variant,
        "gt": str(gt),
        "pred": str(pred),
        "gt_found": int(Path(gt).exists()),
        "pred_found": int(Path(pred).exists()),
    })

for scan in ["EFRGP01_00", "EFRGP01_01"]:
    for crop, (x, y, z) in CROPS.items():
        base, gray, gt, full, raw = crop_paths(scan, crop, x, y, z)

        if not gray.exists():
            raise RuntimeError(f"missing gray: {gray}")
        if not gt.exists():
            raise RuntimeError(f"missing gt: {gt}")
        if not full.exists():
            raise RuntimeError(f"missing full BA-TOS labels: {full}")
        if not raw.exists():
            raise RuntimeError(f"missing raw BA-TOS labels: {raw}")

        add_manifest(scan, crop, "Full BA-TOS", gt, full)
        add_manifest(scan, crop, "w/o size filter", gt, raw)

        nb_raw, nb_filtered = run_variant_no_boundary(scan, crop, gray)
        add_manifest(scan, crop, "w/o boundary rejection", gt, nb_filtered)

        nf_raw, nf_filtered = run_variant_no_flattening(scan, crop, gray)
        add_manifest(scan, crop, "w/o flattening", gt, nf_filtered)

manifest_csv = TABLEDIR / "full_ablation_manifest.csv"
with manifest_csv.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["scan", "crop_id", "variant", "gt", "pred", "gt_found", "pred_found"])
    w.writeheader()
    w.writerows(manifest)

missing = [r for r in manifest if r["gt_found"] == 0 or r["pred_found"] == 0]
print("wrote:", manifest_csv)
print("missing:", len(missing))
for r in missing:
    print(r)

if missing:
    raise SystemExit("Fix missing files before metrics.")
