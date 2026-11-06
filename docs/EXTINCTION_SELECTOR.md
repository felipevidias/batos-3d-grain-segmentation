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
