# Phase 5c — extinction-based marker selection (development crop only)

Date 2026-11-06. Code: `src/batos_v2/extinction_selector.py` (selector), `src/batos_v2/extinction_run.py` (driver),
commit 790b0a6. Frozen code untouched. Crop: **EFRGP01_00_crop04 only** (development crop). Results:
`results/extinction/EFRGP01_00_crop04_summary.csv`, `results/extinction/r_est_vs_gt_radius.csv`, 6 rows in
`results/all_runs.csv` (`phase5c_extinction_*`). Runs: `C:\BA_TOS_DATA\runs\extinction\<variant>\EFRGP01_00_crop04\`.

## 1. What replaces what

Only the `true_dual` stage (frozen weighted selector, `batos_tos_seed_true_dual.py`, K = 3600) is replaced. Everything
else is the frozen pipeline: upstream outputs (openrec r = 15, flatten λ = 1, legacy seed, background, legacy markers)
are read from the frozen run, and downstream e10 spatialisation → e11 ToS-MIN marker builder (NMS 6, EDT weight 0.20)
→ EDT-Otsu control → fixed-mask watershed → E23c gate → selective watershed → `evaluate_v2` (markers, rule
dominant, seed-ignore, no size filter = Table II protocol) are the unmodified frozen calls (the driver reuses
`batos_v2.sensitivity.CropRun`).

## 2. Algorithm (best-effort definitions, documented choices)

1. ToS of the flattened selector image (`hg.component_tree_tree_of_shapes_image3d` on the uint8 image, as in the
   frozen selector and e10; e10 re-checks every node id, altitude and parent altitude and accepted all of them).
2. Dark extrema = Higra extrema (`attribute_extrema`) that are MIN-like (alt < alt(parent)).
3. Boundary rule: reject extrema whose shape touches any of the 6 faces. In practice **0 of 238,168 dark extrema
   touch a face** (Higra's ToS puts voxels on the faces into large border-connected shapes). The rule takes effect in
   step 6, where representative shapes may not grow into a face-touching node.
4. Extinction values computed **directly on the ToS**. Higra's `attribute_extinction_value` / `attribute_dynamics`
   assume monotone altitudes (min/max trees), so they do not apply. Branch rule for an attribute A that is increasing
   on the tree: at each node, among the children that contain a dark extremum, the one with the largest A wins (ties:
   smallest index); the others are branch tops and die there. ext(m) = A(top of m's branch).
   - volume: A = `hg.attribute_volume` (Σ area·|alt − alt(parent)|)
   - dynamics: A = grey-level range (max − min) inside the shape (ToS analogue of height/dynamics)
   - product: ext_volume(m) · ext_dynamics(m); branch membership = winner in both competitions
5. Rank by extinction (ties: darker, then node index), keep the top K = 3600.
6. Representative node of each selected minimum m = the highest MIN-like node on m's own branch whose shape
   contains no other selected minimum and touches no face. These shapes are pairwise disjoint, so no hierarchical
   suppression is needed.
7. CSV `score` = rank score (K − rank + 1)/K, because the frozen e11 builder min-max normalises `score` before mixing it
   with the EDT term and raw extinction values are heavy-tailed. The raw value is kept in the `extinction` column.

## 3. Results — EFRGP01_00_crop04 (272 valid GT grains)

| selector | pool | selected | rep. area med / q75 | ToS markers | grains seeded | split grains (≥ 2 markers) | gate (mv / qv) | P@0.50 | R@0.50 | **F1@0.50** | F1@0.75 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| frozen weighted (Table II) | 148,487 MIN cand. | 1,800 MIN (+1,800 MAX) | 4 / 19 | 398 | 268 | 0 | off (43 / 49) | — | — | **0.985** | 0.978 |
| extinction: volume | 238,168 | 3,600 | 28 / 43 | 477 | 271 | 47 | off (45 / 168) | 0.833 | 0.993 | **0.906** | 0.906 |
| extinction: dynamics | 238,168 | 3,600 | 14 / 28 | 876 | 272 | 192 | off (17 / 123) | 0.441 | 0.993 | **0.611** | 0.604 |
| extinction: product | 238,168 | 3,600 | 18 / 28 | 590 | 268 | 109 | off (21 / 95) | 0.636 | 0.956 | **0.764** | 0.758 |

The gate does not fire, so BA-TOS_sel = BA-TOS_ToS for every variant (identical rows in the summary CSV). "Grains
seeded" and "split grains" use GT for diagnosis only.

**Target ≥ 0.98 F1@0.50: not met by any variant** (best: volume 0.906, −0.079 vs frozen).

### Why: over-seeding, not missed grains
- Recall is as high as or higher than frozen (volume and dynamics seed 271–272 / 272 grains). The loss is entirely in
  precision.
- The e11 builder keeps one marker per internal-dark component, and many grains contain several such components.
  Frozen markers land in exactly one component per seeded grain (**0 split grains**). Extinction ranking also picks
  deep secondary minima inside grains, which seed the extra components and split grains: 47 (volume), 109 (product),
  192 (dynamics). Markers in excluded or face-touching grains are ignored by the seed-ignore protocol (130 frozen,
  153–264 here).
- Dynamics alone ranks contrast and is blind to size, so it promotes small, deep, noisy minima — the worst variant.
  Volume weights by size and is closest to frozen. The frozen score's area prior (`area_preferred_position`, the
  most influential parameter in Phase 5a group 1) is exactly the size control that extinction ranking lacks.

### Memory and runtime
Largest stage peak 3.13 GB (e10 spatialise), selector 3.07 GB, driver 0.62 GB; all stages ran in sequence, so this
work peaked at 3.1 GB (< 6 GB budget). Selector 103–117 s; full downstream ≈ 245–260 s per variant.

## 4. Parameter count of the selection stage

| stage | free parameters |
|---|---|
| frozen weighted selector (true_dual) | **17** swept in Phase 5a (14 influential); `w_support` proven dead |
| extinction selector | **1** numeric (K = 3600, kept equal to frozen) + 1 categorical (variant) |

Structural choices that are fixed (not tuned) but should be stated in the paper: the winner rule (largest-A child), the
representative rule (highest exclusive MIN-like face-free node on the branch), the rank-based score, and the
dynamics definition (grey range of the shape).

## 5. Side check for Phase 5b: e13 r_est against GT radius (GT used for diagnosis only)

| crop | r_est (e13) | openrec r | e13 est. #grains | valid GT | median r_eq valid GT (IQR) | GT / r_est |
|---|---|---|---|---|---|---|
| EFRGP01_00_crop04 | 10.20 | 15 | 1180 | 272 | 14.32 (14.30–14.34) | 1.40 |
| EFRGP01_00_crop05 | 9.00 | 14 | 980 | 129 | 14.35 (14.32–14.39) | 1.59 |
| EFRGP01_01_crop05 | 8.00 | 12 | 1409 | 128 | 14.36 (14.32–14.38) | 1.80 |

r_eq = (3V / 4π)^{1/3}. The GT grains are nearly monodisperse (IQR < 0.1 voxel), and all three crops have the same
true radius (≈ 14.3). Yet e13 underestimates it, and **more on the two crop05 (−37 %, −44 %) than on crop04 (−29 %)**:
on the crops where the gate fires, the GT-free scale estimate is the most wrong. Its derived grain count is also 4–11×
the valid GT count. For Phase 5b: an under-estimated r_est shrinks the openrec radius (14 / 12 vs 15). Group 2 showed
that the openrec radius is the most sensitive single setting (factor 0.75 collapses every crop). A scale-consistency check
(r_est vs grain spacing, or vs the marker count) is a GT-free candidate for the gate.

## 6. Next steps (not run)
- Size-aware extinction: rank by volume extinction but take as representative the extinction node at the scale
  r_est (one extremum per basin of radius ≈ r_est), or merge selected minima whose representatives share an internal
  component before e11.
- K matched to the frozen MIN count (1800) as a GT-free control.
- Only after a variant reaches the target on crop04: the 12 development crops, then frozen before Test.

## 7. Follow-up A (2026-11-06): uniqueness and scale window — target not met

Code commit 37c9b83 (`extinction_selector.py` variants `volume_unique`, `volume_unique_scale`; `extinction_run.py`
multi-crop). Crops: crop04 plus both crop05. Results: `results/extinction/summary.csv` (supersedes
`EFRGP01_00_crop04_summary.csv`, adds grains split), rows in `results/all_runs.csv` (`phase5_extinction_*`).
Both variants are ranked by volume extinction, and both rules are applied **before** the top-K = 3600 cut.

- **volume_unique:** one minimum per 6-connected component of the internal-dark mask (`--internal-source`, as in e10).
  A minimum's component = the component holding most of its extremum flat zone. This is a pre-selection stand-in for
  its exclusive support, which only exists once the selection is made. Minima outside every component are dropped,
  because the e11 builder discards them anyway.
- **volume_unique_scale:** as above, but first keeps only minima whose extinction region has a voxel count in
  [0.05, 2.0] × (4/3)πr³.
  - **Units:** `hg.attribute_volume` is grey-weighted (voxels × grey levels), while (4/3)πr³ is in voxels, so the window is applied
    to the voxel count of the extinction region (the branch top of the volume competition). The literal reading,
    a window on the grey-weighted value, was also run on crop04 (`_greyvol`).
  - **Radius:** r = e13 r_est (as specified), and also r_ps_peak from `docs/SCALE_ESTIMATOR.md`.

### crop04 (272 valid grains; ToS = sel, gate off in every run)

| variant | selected | markers | seeded | split | P@0.50 | R@0.50 | **F1@0.50** | F1@0.75 |
|---|---|---|---|---|---|---|---|---|
| frozen weighted | 1,800 MIN | 398 | 268 | 0 | 0.993 | 0.978 | **0.985** | 0.978 |
| volume (v1) | 3,600 | 477 | 271 | 47 | 0.833 | 0.993 | **0.906** | 0.906 |
| volume_unique | 3,600 of 10,720 components | 1,307 | 272 | 244 | 0.300 | 0.989 | **0.460** | 0.452 |
| volume_unique_scale, r_est 10.20 | 42 (62 in window) | 39 | 36 | 0 | 0.278 | 0.037 | **0.065** | 0.058 |
| volume_unique_scale, r_ps_peak 10 | 60 (83 in window) | 57 | 51 | 0 | 0.353 | 0.066 | **0.111** | 0.093 |
| volume_unique_scale, r_est, literal grey-volume window | 708 (3,394 in window) | 303 | 176 | 30 | 0.587 | 0.449 | **0.508** | 0.450 |

### crop05 (129 / 128 valid grains)

| variant | crop | markers | seeded / split (sel) | gate (mv / qv) | ToS F1@0.50 | sel F1@0.50 | frozen sel |
|---|---|---|---|---|---|---|---|
| volume_unique | 00_crop05 | 1,620 | 129 / 89 | fires (1 / 2) | 0.016 | 0.215 | 0.260 |
| volume_unique | 01_crop05 | 810 | 128 / 69 | fires (1 / 2) | 0.067 | 0.620 | 0.807 |
| volume_unique_scale (r_est and r_ps_peak) | both | — | — | — | **breaks** | **breaks** | |

"Breaks" means **0 candidates**. Only 0–46 minima fall in the scale window and none of them lies in an internal-dark
component, so the frozen e10 fails on the empty CSV (`KeyError: has_reconstructed_region`). This is deterministic,
recorded as `FAILED`, and not re-run.

### Reading
1. **Uniqueness makes it worse, not better (0.906 → 0.460).** The internal-dark mask is fragmented: 19,230
   components on crop04, 17,120 / 15,596 on the crop05 pair, against ~272 / 129 grains. Applying one-per-component
   before the top-K cut spreads the 3,600 candidates over 3,600 different fragments. The e11 builder then emits one marker
   per fragment (1,307 markers, 244 split grains). The frozen selector avoids splits because its candidates concentrate
   in the core fragment of each grain; uniqueness over fragments removes exactly that concentration.
2. **The scale window removes almost everything.** The extinction region of most minima is far smaller than a grain:
   only 62–83 of 238,168 minima on crop04 reach 5 % of a grain volume. In this ToS the grain-core minima lose their
   branch early, against noise siblings inside the grain, so the extinction-region size is not a grain size. Recall
   collapses on crop04, and crop05 has no candidates at all.
3. The literal grey-volume window keeps 3,394 minima but still mixes units. It lands in between (0.508).
4. None of the extinction variants reaches 0.98. The best remains plain volume extinction (0.906). Over-seeding
   (split grains) is the failure mode of every variant with high recall, and under-seeding of every variant that
   avoids splits. The frozen area prior is doing work that these rules do not replace.

### Possible next step (not run)
Make the uniqueness unit the **grain**, not the internal fragment. For example, group fragments by connectivity within
distance r_ps_peak, or merge selected minima whose representative regions lie within one inscribed radius (≈ 10.6)
of each other. Freeze the result on crop04 before touching other crops.
