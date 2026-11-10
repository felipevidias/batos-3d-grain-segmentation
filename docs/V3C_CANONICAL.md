# BA-TOS v3c (canonical-scale normalization) — design record, Phase 5i

Written 2026-11-10 **before any v3c run**. Design data: the 12 EFR dev crops (f = 1.0) and the synthetic scale set
(f = 0.7 / 1.4 / 2.0, `data/metadata/crops_synth_scale.json`). Test set 1 and test set 2 are not read during design.
Branch `journal-extension-v32` (worktree `..\batos-v32`).

## Idea

v3 rewrites every v2 constant in grain units (s = r_granulo / 10). v3c instead resamples the **image**, so its grains
have the canonical radius R_REF = 10 voxels (the dev scale). It then runs the frozen v2-fixflat pipeline with every
constant at its frozen value, and maps the result back to the native grid. r_granulo is the only scale input.

## Pipeline (`src/batos_v3c/canonical.py`; GT never read)

1. **Scale.** r_granulo = `batos_v3.scale.estimate(gray)` on the native crop. This is the same estimator as v3 (integer,
   2..40). f = R_REF / r_granulo.
2. **Resample in.** If f = 1, the native file is used unchanged (no interpolation). Otherwise:
   - grey = `scipy.ndimage.zoom(gray, Nc/N, order=3, grid_mode=True, mode="grid-mirror")` (cubic, voxel-edge
     aligned, per-axis factor Nc/N with Nc = round(N f));
   - the result is clipped to the input dtype range, rounded, and cast back to the input dtype.
3. **Canonical run.** `batos_v3.pipeline.run_crop(..., preset="fixflat")`, the code path that reproduces
   batos-paper-v2-fixflat bit-for-bit. The e13 scale is replaced by a params file with estimated_grain_radius = 10 and
   openrec_radius = max(8, round(1.5 · 10)) = 15: the frozen rule evaluated at the canonical radius. All other flags are
   the frozen `BASE_FLAGS` plus flatten `--variant A`. So the Otsu-mask min component (500), the closing (1), K = 3600,
   NMS 6 and the gate (3 / 4) all act at canonical scale.
4. **Size filter "own".** 500 voxels, applied to the canonical label volume before mapping back. This is equivalent
   to 500 / f³ native voxels.
5. **Map back.**
   - Labels: nearest neighbour. Native voxel i takes canonical voxel floor((i + 0.5) · Nc/N).
   - Markers (needed for seed-ignore evaluation): each canonical marker voxel j maps to native
     round((j + 0.5) · N/Nc − 0.5), keeping its label. The prediction ↔ seed link is by label value, as in e09b.
6. **Failure.** A stage failure is recorded as FAILED.json and scored as an empty segmentation (same as 5e/5f).

Identity (checked before the dev run): on a crop with r_granulo = 10, v3c is fixflat with openrec radius 15. On
EFRGP01_00_crop04 (frozen e13 radius 15) it must reproduce the existing batos-paper-v2-fixflat segmentation
bit-for-bit.

Note: on dev, v3c differs from fixflat only on the crops whose e13 radius ≠ 15 (e13 r̂ 8–9 → 12/14; 10.39 → 16).

## Evaluation

Same protocol as Phase 5e step 3 (`run_v3_scale_test.evaluate_crop`):
- native GT; valid GT ≥ round(1000 f³) voxels and not touching a face;
- routes ToS / sel; size filter none / own;
- modes markers_dominant (seed-ignore, primary) and no_markers_fraction050;
- IoU 0.25 / 0.50 / 0.75.

v3 and fixflat are compared from their existing 5e evaluations.

Also reported per factor:
- markers per grain (sel markers in native space; all GT and valid GT), seeded and split fractions;
- boundary cost of the resampling against v3, on valid GT grains matched (IoU > 0.5) by **both** v3c and v3: mean
  matched IoU and median signed / absolute relative volume error.

## Acceptance (fixed before running)

Primary metric: crop-macro F1@0.50, BA-TOS_sel, seed-ignore, no size filter.

1. f = 1.0: ≥ 0.976 (v2-fixflat dev).
2. "≈ flat" on the synthetic factors, made concrete here: f = 0.7, 1.4 and 2.0 each ≥ 0.966 (the f = 1.0 floor
   − 0.010). For reference, v3 scores 0.975 / 0.995 / 0.938, so v3 would fail this at f = 2.0.

Both must hold. If they do:
- freeze as `batos-paper-v3c` (config, tag, SHA-256 of executed code);
- cut test set 2 (seed 20261109, same GT-free rule as test set 1, 52 crops, no overlap with dev or test set 1);
- run v3 and v3c once;
- report per specimen (EOR/EOS emphasised) and the boundary cost.

If not: report and stop. Test set 2 stays unseen.
