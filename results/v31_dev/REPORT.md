# Phase 5f step 2 — v3.1 validation (dev + synthetic scale set): FAILED, not frozen

v3.1 = v3 with count constants (kind N) set to v2 · max(1, s⁻³), floored at the v2 capacity instead of v2/s³.
N constants: `tos_seed --max-selected-nodes` 900, `markers_legacy --max-markers` 420, `true_dual --max-selected-nodes`
3600 (one flag covers K = 1800 MIN + 1800 MAX). At s = 2 the values are 900 / 420 / 3600, against 112 / 52 / 450 in v3.
Design data: the 12 EFR dev crops plus the synthetic f = 0.7/1.4/2.0 set. The test crops were not touched.

Full tables: `tables.md`. CSVs: `crop_macro.csv`, `pooled.csv`, `per_crop.csv`, `seed_density.csv`.

## Acceptance criteria (Phase 5f step 2)

| criterion | v3 | v3.1 | result |
|---|---|---|---|
| dev f = 1.0 stays 0.989 | 0.989 | 0.989 (identical flags on all 24 crops with s ≤ 1; crop04 re-run bit-identical, `identity_check.json`) | **pass** |
| f = 1.4 holds or improves | 0.995 | 0.941 (2 crops < 0.90: both crop05, 0.739 / 0.735) | **fail** |
| f = 2.0 holds or improves | 0.938 | 0.326 (12/12 crops < 0.90) | **fail** |
| the two crop05 at f = 2.0 (0.667 / 0.621) improve | 0.667 / 0.621 | 0.247 / 0.290 | **fail** |

Crop-macro F1@0.50, BA-TOS_sel, seed-ignore, no size filter. Pooled: 0.995 → 0.951 (f = 1.4), 0.947 → 0.315 (f = 2.0).

## Seed density per GT grain (medians; GT counts used for diagnosis only)

| f | v3 cand./grain | v3.1 cand./grain | v3 markers/grain | v3.1 markers/grain | v3 pred/grain | v3.1 pred/grain |
|---|---|---|---|---|---|---|
| 0.7 | 3.87 | 3.87 | 0.75 | 0.75 | 0.75 | 0.75 |
| 1.0 | 3.02 | 3.02 | 0.62 | 0.62 | 0.62 | 0.62 |
| 1.4 | 2.40 | 6.58 | 0.52 | 0.53 | 0.52 | 0.53 |
| 2.0 | 1.95 | 15.58 | 0.42 | 2.62 | 0.42 | 2.62 |

## Why it fails

Holding capacity at the v2 value when grains are larger gives about 15 MIN candidates per grain at f = 2.0. The
downstream marker refinement no longer reduces those to one marker per grain: there are 2.6 markers per grain
(0.42 in v3), and each extra marker seeds its own watershed basin.

The failure is pure over-segmentation into small fragments:
- FN = 0 on every crop at f = 1.4 and 2.0, and recall is 1.00 everywhere.
- At f = 2.0 there are on average 247 predictions per crop and 68 FP@0.50. With the crop's own size filter
  (500 · s³ voxels) this drops to 51 predictions and 0.17 FP, and F1 rises to 0.996 (v3: 1.000).
- So the extra markers split grains into sub-filter fragments. Precision without the filter is 0.14–0.32.

At f = 1.4 the effect appears only on the two crop05. Those are the gate crops with the weakest support, at 0.86
markers per grain against about 0.52 elsewhere. The gate does not fire on any v3.1 crop at f = 2.0 (region-voxel
median 47–57 against limits 24–28), so the EDT-Otsu fallback does not catch the fragmentation.

The capacity cap in v3 was doing useful work. It is the de facto limit that keeps seed multiplicity down at large
scale. The test-set-1 EO failure (too few candidates) and this dev failure (too many markers) pull in opposite
directions, so a pure capacity rule cannot satisfy both without a per-grain merge or suppression step.

## Consequence

Steps 3 and 4 were **not** executed: no freeze, no tag `batos-paper-v3.1`, no test set 2 defined or cut. Test set 2
(seed 20261109) remains unseen. The step-3/4 scripts (`freeze_v31.py`, `define_test_crops2.py`, `cut_test2.py`,
`run_test2.py`, `aggregate_test2.py`) are committed but have not run. They were dry-tested on test-1 data only.
Test set 1 numbers for v3 are unchanged (`results/test_v3/REPORT.md`).
