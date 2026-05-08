#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
PY="${PY:-/home/felipe/Pesquisa-Grain_Seg/Filtro_Conexo/venv/bin/python3}"

GRAY="/mnt/windows_p3/Users/USER/Downloads/PhDUpload/PhDUpload/EFR/Greyscale/EFRGP01_01.tif"
LABEL="/mnt/windows_p3/Users/USER/Downloads/PhDUpload/PhDUpload/EFR/Label/extracted_EFRGP01_01/EFRGP01_01.tif"

OUTROOT="$ROOT/results_batos/EFRGP01_01_selected_crops"
mkdir -p "$OUTROOT"

SUMMARY="$OUTROOT/EFRGP01_01_selected_crops_methods_summary.csv"
echo "scan,crop_id,method,gt_valid,pred_labels,pred_eval,pred_ignored,precision050,recall050,f1050,fp050,fn050,overseg_ratio" > "$SUMMARY"

append_instance_summary() {
  local CROP="$1"
  local METHOD="$2"
  local CSV="$3"

  "$PY" - <<PY
import csv
from pathlib import Path

csv_path = Path("$CSV")
with open(csv_path) as f:
    rows = list(csv.DictReader(f))

if not rows:
    raise SystemExit(f"Empty metrics file: {csv_path}")

r = rows[0]

line = [
    "EFRGP01_01",
    "$CROP",
    "$METHOD",
    r["gt_valid"],
    r["pred_total"],
    r["pred_eval"],
    r["pred_ignored_border"],
    r["precision"],
    r["recall"],
    r["f1"],
    r["fp"],
    r["fn"],
    r["overseg_ratio_pred_eval_gt"],
]

with open("$SUMMARY", "a", newline="") as f:
    csv.writer(f).writerow(line)
PY
}

run_one() {
  local CROP="$1"
  local X="$2"
  local Y="$3"
  local Z="$4"
  local SIZE="$5"

  local CROP_DIR="$OUTROOT/${CROP}_x${X}_y${Y}_z${Z}_s${SIZE}"
  local INPUT_DIR="$CROP_DIR/00_input"

  mkdir -p "$INPUT_DIR"

  local PREFIX="EFRGP01_01_${CROP}_x${X}_y${Y}_z${Z}_s${SIZE}"
  local GRAY_CROP="$INPUT_DIR/${PREFIX}_gray_crop_8u.tif"
  local LABEL_CROP="$INPUT_DIR/${PREFIX}_label_crop.tif"

  echo ""
  echo "===================================================="
  echo "[EFRGP01_01] $CROP | x=$X y=$Y z=$Z size=$SIZE"
  echo "===================================================="

  if [ ! -f "$GRAY_CROP" ] || [ ! -f "$LABEL_CROP" ]; then
    "$PY" "$ROOT/code/segmentation/batos/batos_make_paired_crop.py" \
      --gray "$GRAY" \
      --label "$LABEL" \
      --out-dir "$INPUT_DIR" \
      --x "$X" \
      --y "$Y" \
      --z "$Z" \
      --size "$SIZE" \
      --prefix "$PREFIX"
  else
    echo "[crop] already exists: $GRAY_CROP"
  fi

  echo ""
  echo "--------------------"
  echo "[1/3] BA-TOS"
  echo "--------------------"

  "$ROOT/code/segmentation/batos/run_batos_pipeline.sh" \
    "$GRAY_CROP" \
    "$CROP_DIR/batos_pipeline" \
    15 \
    500

  "$PY" "$ROOT/code/segmentation/batos/batos_instance_metrics.py" \
    --gt "$LABEL_CROP" \
    --pred "$CROP_DIR/batos_pipeline/05_watershed/batos_distance_filtered_min500_labels.tif" \
    --out-dir "$CROP_DIR/batos_eval_instance" \
    --min-gt-area 1000 \
    --ignore-gt-border \
    --thresholds 0.50 \
    --name "${CROP}_batos"

  append_instance_summary \
    "$CROP" \
    "BA-TOS" \
    "$CROP_DIR/batos_eval_instance/${CROP}_batos_instance_metrics.csv"

  echo ""
  echo "--------------------"
  echo "[2/3] Otsu+WS"
  echo "--------------------"

  "$PY" "$ROOT/code/segmentation/batos/baseline_otsu_distance_watershed.py" \
    --image "$GRAY_CROP" \
    --out-dir "$CROP_DIR/otsu_ws" \
    --prefix "${CROP}_otsu_ws" \
    --threshold-offset -10 \
    --min-distance 8 \
    --close-radius 1

  "$PY" "$ROOT/code/segmentation/batos/batos_instance_metrics.py" \
    --gt "$LABEL_CROP" \
    --pred "$CROP_DIR/otsu_ws/${CROP}_otsu_ws_labels.tif" \
    --out-dir "$CROP_DIR/otsu_eval_instance" \
    --min-gt-area 1000 \
    --ignore-gt-border \
    --thresholds 0.50 \
    --name "${CROP}_otsu"

  append_instance_summary \
    "$CROP" \
    "Otsu+WS" \
    "$CROP_DIR/otsu_eval_instance/${CROP}_otsu_instance_metrics.csv"

  echo ""
  echo "--------------------"
  echo "[3/3] Macke-MinTree+WS"
  echo "--------------------"

  "$ROOT/code/segmentation/batos/run_macke_mintree_one_crop.sh" \
    "$GRAY_CROP" \
    "$LABEL_CROP" \
    "$CROP_DIR/macke_mintree"

  append_instance_summary \
    "$CROP" \
    "Macke-MinTree+WS" \
    "$CROP_DIR/macke_mintree/04_instance/macke_mintree_instance_metrics.csv"
}

run_one "crop02" 600 800 300 200
run_one "crop03" 900 600 300 200
run_one "crop04" 500 500 700 200
run_one "crop05" 800 900 700 200
run_one "crop06" 350 650 500 200
run_one "crop07" 700 400 500 200

echo ""
echo "===================================================="
echo "[EFRGP01_01] Final summary"
echo "===================================================="
cat "$SUMMARY"
