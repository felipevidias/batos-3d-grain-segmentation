# Phase 5b side study — GT-free scale estimator and the crop05 failure

Date 2026-11-06. Crops: EFRGP01_00_crop04, EFRGP01_00_crop05, EFRGP01_01_crop05 only. Code (commit 37c9b83):
`src/batos_v2/scale_granulo.py` (estimators), `src/batos_v2/openrec_lowmem.py` (exact low-memory openrec),
`src/batos_v2/scale_rerun.py` (frozen pipeline with a new openrec radius). Tables: `results/scale/scale_estimates.csv`,
`results/scale/<crop>_granulo.json` (full pattern spectra), `results/scale/rerun_summary.csv`; 6 rows in
`results/all_runs.csv` (`phase5_frozen_openrec_from_r_ps_peak/*`). Runs: `C:\BA_TOS_DATA\runs\scale_granulo\r_ps_peak\`.
Frozen code untouched; every stage peaked ≤ 3.66 GB (sequential, < 6 GB budget).

## 1. Estimators (GT never read)

Foreground = e13's own Otsu foreground (robust p1–p99 uint8 → Otsu ≥ → closing ball(1) → fill holes), imported read-only from
the frozen `e13_estimate_scale.py`.

- **Pattern spectrum.** Binary openings by balls of radius r = 2…31, computed exactly with EDTs (erosion = {EDT > r},
  dilation of that core = {distance ≤ r}). The crop outside is not treated as background. PS(r) = V(r) − V(r+1). Two
  estimates: `r_ps_peak` = argmax PS, and `r_ps_median` = radius at which half the removed volume is gone.
- **EDT local maxima.** Median EDT over the 3×3×3 local maxima inside the foreground, with no percentile cut and no top-N
  truncation. (e13 itself is a variant of this: 5³ window, peaks ≥ EDT p75, the 1,500 largest.)

## 2. Estimates against GT (GT used for diagnosis only)

| crop | e13 r_est | **r_ps_peak** | r_ps_median | r_EDT-localmax | Otsu fg | GT label fraction | valid GT | GT median r_eq | **GT median inscribed r** |
|---|---|---|---|---|---|---|---|---|---|
| 00_crop04 | 10.20 | **10** | 10 | 2.00 | 0.655 | 0.675 | 272 | 14.32 | **10.63** |
| 00_crop05 | 9.00 | **11** | 10 | 2.24 | 0.374 | 0.324 | 129 | 14.35 | **10.63** |
| 01_crop05 | 8.00 | **11** | 10 | 2.24 | 0.378 | 0.328 | 128 | 14.36 | **10.63** |

(r_eq = (3V/4π)^{1/3}; inscribed r = maximum EDT inside each valid GT grain.)

**Correction to the earlier side check.** Earlier I compared r_est with the GT *equivalent* radius (14.3). The grains are
not spherical: their median **inscribed** radius is **10.63 on all three crops** (r_eq / r_in ≈ 1.35). A ball-based
estimator (EDT, openings) measures the inscribed radius, so that is the right yardstick:
- e13 is about right on crop04 (10.20, −4 %) but underestimates on both crop05: **9.0 (−15 %) and 8.0 (−25 %)**. The cause
  is the crop, not the grains. Only 37 % of crop05 is foreground (66 % on crop04), and 16 % of that foreground carries
  no GT label (bright non-grain material). This changes which EDT peaks pass e13's p75 cut and its top-1,500 truncation.
- The pattern-spectrum peak is stable (10 / 11 / 11, within ±0.6 of 10.63). The spectrum is sharp and vanishes at
  r = 12–13 on every crop, as expected for near-monodisperse grains.
- The EDT local-maxima median is useless (2.0–2.2): without e13's cuts it is dominated by noise maxima.
- Otsu foreground and GT agree: 96–97 % of GT voxels are foreground, and the GT r_eq restricted to the foreground is unchanged.

## 3. Frozen pipeline with openrec radius = max(8, round(1.5 · r_ps_peak))

Only the radius changes (15 / 16 / 16 instead of 15 / 14 / 12). The openrec stage runs `batos_v2.openrec_lowmem`,
an exact re-implementation: a ball erosion split into 2-D disc erosions per z-shift, with scipy's reflect boundary
applied axis by axis. It is bit-identical to `skimage.erosion(img, ball(r))` on random volumes (r = 5, 9, 12), and the
whole stage is bit-identical to the frozen openrec on crop04 at r = 15 (4/4 outputs). Its peak is 740 MB instead of
~4.8 GB at r = 15; the frozen scipy erosion needs ~18 GB at r = 20. All stages, including crop04, were re-run from openrec.

| crop | route | openrec r | gate (mv / qv) | markers | grains seeded / split | F1@0.50 | F1@0.75 | frozen F1@0.50 (markers) |
|---|---|---|---|---|---|---|---|---|
| 00_crop04 | ToS = sel | 15 | off (43 / 49) | 398 | 268 / 0 | **0.985** | 0.978 | 0.985 (398) — identical |
| 00_crop05 | ToS | 16 | fires (2 / 4) | 254 | 90 / 41 | 0.007 | 0.007 | 0.012 (53) |
| 00_crop05 | sel | 16 | fires | 254 | 129 / 55 | **0.722** | 0.462 | 0.260 (53) |
| 01_crop05 | ToS | 16 | fires (2 / 3) | 364 | 104 / 66 | 0.016 | 0.016 | 0.066 (191) |
| 01_crop05 | sel | 16 | fires | 364 | 128 / 55 | **0.702** | 0.511 | 0.807 (191) |

**crop05 does not recover. The scale estimate is not the root cause.**
- With the correct radius the ToS support still collapses (mv 2, qv 3–4), and the gate still fires on both crop05. ToS
  F1 stays at ~0.
- What changes is the **number** of ToS markers (53 → 254, 191 → 364). The EDT-Otsu fallback is count-matched to that
  number (`--match-count-from`), so the fallback gains seeds:
  - 00_crop05 was starved (48 grains seeded) and rises to 0.722.
  - 01_crop05 becomes over-seeded (55 split grains) and drops to 0.702.
  This matches Phase 5a group 2 (`openrec_radius_floor = 16`: 00_crop05 sel 0.26 → 0.70).
- The collapse lives in the selector or markers on these crops, not in the openrec scale. A candidate cause is the
  composition of crop05 (37 % foreground, unlabelled bright material). Phase 5b should test it, with the count-matching
  of the fallback as a separate factor.

## 4. For Phase 5b
- Replace e13 by `r_ps_peak` as the scale estimate. It is more stable across crops and needs no cut or top-N
  parameter; its single remaining choice is the radius range 2–30.
- Decouple the fallback marker count from the collapsed ToS count. For example, use count = foreground volume /
  (4/3 π r³ · shape factor), or let EDT-Otsu choose its own count. On 00_crop05 the fallback is limited by the count,
  not by the method.
