#!/usr/bin/env bash
set -euo pipefail

# BA-TOS: Boundary-Aware ToSConOp-based Marker Extraction
#
# Final pipeline:
# 1. Conservative flattening
# 2. ToSConOp/autodual candidate extraction
# 3. Boundary-connected dark background rejection
# 4. Internal marker refinement
# 5. BA-TOS distance watershed
#
# Uso:
# ./run_batos_pipeline.sh INPUT_GRAY.tif OUTDIR [LAMBDA]

INPUT="${1:?uso: INPUT_GRAY.tif OUTDIR [LAMBDA]}"
OUTDIR="${2:?uso: INPUT_GRAY.tif OUTDIR [LAMBDA]}"
LAMBDA="${3:-15}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_ABS="$(realpath "$INPUT")"
OUT_ABS="$(realpath -m "$OUTDIR")"

STEM="$(basename "$INPUT_ABS")"
STEM="${STEM%.tif}"

mkdir -p "$OUT_ABS"/{00_input,01_preprocess,02_tos_seed,03_background,04_markers,05_watershed,06_qc}

echo ""
echo "======================================"
echo "[BA-TOS] input  = $INPUT_ABS"
echo "[BA-TOS] outdir = $OUT_ABS"
echo "[BA-TOS] lambda = $LAMBDA"\necho "[BA-TOS] min label size = ${4:-500}"
echo "======================================"

echo ""
echo "======================================"
echo "[1/5] Conservative flattening"
echo "======================================"

python3 "$SCRIPT_DIR/batos_preprocess.py" "$INPUT_ABS" \
  --out-dir "$OUT_ABS/01_preprocess" \
  --prefix "$STEM" \
  --lambda-xy "$LAMBDA"

FLAT="$(find "$OUT_ABS/01_preprocess" -maxdepth 1 -type f -name "${STEM}_flatten_lambda*_8u.tif" | sort | tail -n 1)"

if [ ! -f "$FLAT" ]; then
  echo "[BA-TOS][ERRO] Flatten não encontrado."
  exit 1
fi

echo "[BA-TOS] flattened = $FLAT"

echo ""
echo "======================================"
echo "[2/5] Raw ToSConOp/autodual candidate extraction"
echo "======================================"

python3 "$SCRIPT_DIR/batos_tos_seed.py" "$FLAT" \
  --out-dir "$OUT_ABS/02_tos_seed" \
  --height-factor 0.012 \
  --area-factor 22 \
  --area-min 2 \
  --area-max 2200 \
  --delta-parent-min 0 \
  --dark-percentile-max 78 \
  --min-score 0.35 \
  --max-selected-nodes 900 \
  --support-area-min 20 \
  --support-area-max 45000 \
  --support-max-steps 16 \
  --component-min 2 \
  --component-max 1500

RAW_SEED="$(find "$OUT_ABS/02_tos_seed" -maxdepth 1 -type f -name "*_tosconop_autodual_seed_raw.tif" | sort | tail -n 1)"

if [ ! -f "$RAW_SEED" ]; then
  echo "[BA-TOS][ERRO] Seed bruta ToSConOp não encontrada."
  exit 1
fi

echo "[BA-TOS] raw ToS seed = $RAW_SEED"

echo ""
echo "======================================"
echo "[3/5] Boundary-connected dark background rejection"
echo "======================================"

python3 "$SCRIPT_DIR/batos_background.py" "$INPUT_ABS" \
  --flattened "$FLAT" \
  --seed "$RAW_SEED" \
  --out-dir "$OUT_ABS/03_background" \
  --prefix "$STEM" \
  --dark-percentile 40 \
  --otsu-offset -10 \
  --dark-threshold-max 135 \
  --dark-connectivity 1 \
  --seed-area-min 1 \
  --seed-area-max 3000 \
  --max-external-overlap 0.10 \
  --min-internal-overlap 0.00 \
  --internal-search-radius 6 \
  --center-dark-max 210 \
  --min-ring-contrast -5 \
  --max-ring-external-coverage 0.70 \
  --nms-min-dist 5

INTERNAL="$(find "$OUT_ABS/03_background" -maxdepth 1 -type f -name "*_internal_dark_candidates.tif" | sort | tail -n 1)"

if [ ! -f "$INTERNAL" ]; then
  echo "[BA-TOS][ERRO] internal_dark_candidates não encontrado."
  exit 1
fi

echo "[BA-TOS] internal dark candidates = $INTERNAL"

echo ""
echo "======================================"
echo "[4/5] Internal marker refinement"
echo "======================================"

python3 "$SCRIPT_DIR/batos_markers.py" "$INPUT_ABS" \
  --flattened "$FLAT" \
  --internal-dark "$INTERNAL" \
  --out-dir "$OUT_ABS/04_markers" \
  --prefix "$STEM" \
  --connectivity 1 \
  --area-min 3 \
  --area-max 900 \
  --max-dim 30 \
  --max-elongation 8 \
  --min-fill 0.005 \
  --center-dark-max 205 \
  --min-ring-contrast -2 \
  --nms-min-dist 7 \
  --max-markers 420

MARKER_COMPONENTS="$(find "$OUT_ABS/04_markers" -maxdepth 1 -type f -name "*_batos_marker_components_raw.tif" | sort | tail -n 1)"

if [ ! -f "$MARKER_COMPONENTS" ]; then
  echo "[BA-TOS][ERRO] componentes de marker não encontrados."
  exit 1
fi

echo "[BA-TOS] marker components = $MARKER_COMPONENTS"

echo ""
echo "======================================"
echo "[5/6] BA-TOS distance watershed"
echo "======================================"

python3 "$SCRIPT_DIR/batos_watershed_distance.py" \
  --image "$INPUT_ABS" \
  --markers "$MARKER_COMPONENTS" \
  --out-dir "$OUT_ABS/05_watershed" \
  --prefix batos_distance \
  --threshold-offset -10 \
  --close-radius 1 \
  --force-markers-foreground \
  --force-radius 1

RAW_RESULT_LABELS="$OUT_ABS/05_watershed/batos_distance_labels.tif"
RAW_RESULT_COLOR="$OUT_ABS/05_watershed/batos_distance_colored_labels_rgb.tif"
RAW_RESULT_OVERLAY="$OUT_ABS/05_watershed/batos_distance_overlay_rgb.tif"

echo ""
echo "======================================"
echo "[6/6] Size-aware predicted-instance filtering"
echo "======================================"

MIN_LABEL_SIZE="${4:-500}"

python3 "$SCRIPT_DIR/filter_predicted_labels_by_size.py" \
  --labels "$RAW_RESULT_LABELS" \
  --out-dir "$OUT_ABS/05_watershed" \
  --prefix batos_distance_filtered \
  --min-size "$MIN_LABEL_SIZE" \
  --relabel

RESULT_LABELS="$OUT_ABS/05_watershed/batos_distance_filtered_min${MIN_LABEL_SIZE}_labels.tif"
RESULT_COLOR="$OUT_ABS/05_watershed/batos_distance_filtered_min${MIN_LABEL_SIZE}_colored_labels_rgb.tif"
RESULT_OVERLAY="$RAW_RESULT_OVERLAY"

echo ""
echo "======================================"
echo "[BA-TOS] Counts"
echo "======================================"

python3 "$SCRIPT_DIR/count_tiff_labels.py" "$MARKER_COMPONENTS" --binary \
  | tee "$OUT_ABS/06_qc/count_markers.txt"

python3 "$SCRIPT_DIR/count_tiff_labels.py" "$RESULT_LABELS" \
  | tee "$OUT_ABS/06_qc/count_labels.txt"

echo ""
echo "======================================"
echo "[BA-TOS] Open QC"
echo "======================================"
echo "gmic $OUT_ABS/03_background/${STEM}_external_internal_overlay_rgb.tif a z"
echo "gmic $OUT_ABS/04_markers/${STEM}_batos_marker_components_overlay_rgb.tif a z"
echo "gmic $OUT_ABS/04_markers/${STEM}_batos_marker_centroids_overlay_rgb.tif a z"
echo "gmic $RESULT_OVERLAY a z"
echo "gmic $RESULT_COLOR a z"

echo ""
echo "[BA-TOS] Done."
