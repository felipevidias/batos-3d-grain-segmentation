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

## Amendment 1 (2026-11-10, user protocol clarification; before any v3c result)

The canonical **field of view** is held at 200³, so the synthetic comparison isolates grain scale from crop extent.

**Synthetic windows** (`scripts/make_synthetic_canon.py`, `data/metadata/crops_synth_canon.json`):
- The probe is the existing 200³ synthetic crop at the same centre. Its r_granulo is recomputed and must equal the
  5e record.
- Native side S = round(200 · r_probe / 10), re-cut from the full volume with the same transform as
  `make_synthetic_scale.py`. The S = 200 case is checked to reproduce the 5e crop exactly.
- Native sides: f = 0.7 → 140 (r 7); 1.4 → 280 (r 14); 2.0 → 400 (r 20, 10 crops) or 420 (r 21, 2 crops).
  - The user's note gave 143 for 0.7. The stated rule round(200 · r / 10) with r = 7 gives 140, which maps to exactly
    200 at f = 10/7, so 140 is used.
- Every canonical crop is therefore exactly 200³.

**Scale input.** v3c uses r_probe as its scale input (one estimate per window). This keeps the canonical side exactly
200 instead of re-estimating on the larger window. For the dev crops the crop is its own probe.

**Valid GT.** ≥ round(1000 · (r_probe / 10)³) native voxels, i.e. 1000 at canonical scale (8 000 at r = 20, 9 261 at
r = 21, 343 at r = 7, 2 744 at r = 14). Not touching a face. Evaluation stays in native space against native GT.

**References for the synthetic factors.** v3 and fixflat have no runs on these windows. Their 5e numbers (200³ native
windows) are shown for reference only, and the boundary cost against v3 is computed only at f = 1.0 (same window). The
acceptance thresholds are unchanged.

**Test set 2 (if v3c passes).**
- Per test-2 crop, r_granulo is estimated on the 200³ probe.
- v3c then runs on a re-cut native window of side round(200 · r / 10) at the same centre. The sizes used are reported.
- v3 runs on the 200³ probe (its own protocol).
- Both are evaluated against their own window's GT.

## Amendment 2 (2026-11-10, implementation details; written while the dev/synthetic run is in progress)

The acceptance criteria and thresholds above are **unchanged**. Results seen when this was written: the 12 f = 1.0
crops and 5 of the 12 f = 2.0 crops (per-crop log lines only, no aggregate). Nothing below depends on them.

**Evaluator (memory).** On a 400³ window, the frozen `batos_v2.evaluate_v2.evaluate` peaks at 4.8 GB; the driver
reached 5.8 GB on EFRGP01_00_crop02_x2_canon. This is above the 4 GB working rule, so the run was stopped there.
`src/batos_v3c/eval_lean.py` accumulates the same voxel counts over z-slabs and then applies the frozen matching and
metric code (the frozen module is imported, never edited). `scripts/check_eval_lean.py` found 0 differences against
the frozen evaluator:
- 96 metric rows on the 12 dev crops (2 routes × 2 filters × 2 modes), plus seed statistics, boundary matching and an
  empty prediction;
- the 8 rows of the 400³ crop scored with the frozen evaluator (read with round-trip float parsing).

Peak on 400³: 2.1 GB. All later crops are scored with `eval_lean`. The crop02_x2 rows are the frozen evaluator's, and
they are identical.

**Freeze hashes.** `configs/batos_paper_v3_scalefree.yaml` recorded LF bytes, but this worktree checks out CRLF
(core.autocrlf). Seven of its files therefore fail a raw-byte SHA check, even though `git diff` against the tag
`batos-paper-v3-scalefree` is empty and their LF-normalised SHA equals the record. The test-2 driver checks the v3
freeze on LF-normalised bytes **and** requires an empty `git diff` against that tag. The v3c freeze
(`scripts/freeze_v3c.py`) records LF-normalised hashes.

**Test set 2 (only if both criteria hold).**
- Definition: `define_test_crops2.py`, unchanged except its description. Seed 20261109, 52 crops, no overlap with dev
  or test set 1.
- Cut: `cut_test2_v3c.py`. It produces the 200³ probes (unchanged cutter) and, per crop:
  - r_probe = `batos_v3.scale.estimate(probe)`;
  - S = round(200 · r_probe / 10);
  - window origin = probe origin + ⌊(200 − S) / 2⌋ per axis (same centre; exact for even S, which holds for every
    integer r);
  - if the window would leave the volume, it is shifted inward by the minimum amount and the shift is recorded;
  - S larger than the volume fails loudly.
- Run: `run_test2_v3c.py`. v3 runs on the probe (unchanged `run_test_v3.run_v3`); v3c runs on the window with
  r_probe. Each method is evaluated against its own window's GT with `eval_lean`, under both valid-GT rules:
  - `fixed1000`: ≥ 1000 native voxels. Primary for v3, as on test set 1.
  - `canon1000`: ≥ round(1000 (r/10)³) native voxels. Primary for v3c, per Amendment 1.
- Boundary cost on test 2: on GT grains valid (canon1000) in **both** windows and matched with IoU > 0.5 by both
  methods (GT label ids are global to the scan), report the per-crop mean IoU and the median signed and absolute
  relative volume error.
- Aggregate: `aggregate_test2_v3c.py`. Per specimen, EFR_internal, external and all. Crop-macro and pooled tables are
  kept separate. Native sizes are reported. Rows go to results/all_runs.csv.
