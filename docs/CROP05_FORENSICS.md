# Crop05 forensics — why the support collapses on EFRGP01_00_crop05 and EFRGP01_01_crop05

Date 2026-11-06. Read-only analysis of the frozen outputs (`C:\BA_TOS_DATA\runs\batos_paper_v2_frozen`), with
crop04 (EFRGP01_00_crop04) as the reference. Nothing was re-run except two read-only re-executions of frozen functions
(e13 foreground and the λ = 1 flatten), both bit-identical to the frozen outputs. GT is used for diagnosis only. Renders:
`C:\BA_TOS_DATA\runs\crop05_forensics\{EFRGP01_00_crop05,EFRGP01_01_crop05,EFRGP01_00_crop04}_midslices.png`.
Stage tables: `results/crop05_forensics/stage_comparison.csv`.

## Summary

**First divergence: the λ = 1 flatten stage that produces the selector image** (`src/batos/batos_preprocess.py`,
`flatten_volume_u8`). It adds back the offset `bg_ref = median(blur[blur > 0])`, where `blur` is a per-slice Gaussian
(σ = 1). The `> 0` guard is meant to skip empty background, but a float Gaussian leaves tiny positive tails over zero
plateaux. On crop05 the openrec image is 58–67 % exactly zero, from the exterior air plus the pores. As a result, 51 % of
the "positive" blurred voxels are below one grey level and **bg_ref collapses to 0.9, against 186 on crop04**.

After the robust rescale, the selector has lost its offset:
- grain bodies drop to 44 / 33 (crop04: 165);
- grain cores, pores and exterior all sit at 0–7.

Every downstream stage inherits this compressed image:
- the ToS MIN pool shrinks;
- the MIN height scale is 3–4 (crop04: 32);
- the dark prior `dark_alt_max` = 46 / 42 (crop04: 178);
- the internal-dark mask becomes fragments of dim grain bodies (75–83 per grain; crop04: 13);
- 97 % of the selected MIN points land outside the grains;
- only 53 / 191 markers are made, with 2-voxel supports, so the gate fires.

The precondition is geometric: **crop05 straddles the specimen boundary**. About 45 % of the crop is exterior air
beyond the confining membrane.

## (a) What is the 16 % unlabelled bright foreground?

It is the **specimen's confining membrane / sleeve**, not grain material and not a rim around every grain.

| | 00_crop05 | 01_crop05 | 00_crop04 (ref.) |
|---|---|---|---|
| Otsu foreground (e13) | 0.374 | 0.378 | 0.655 |
| GT-labelled fraction | 0.324 | 0.328 | 0.675 |
| unlabelled foreground, share of foreground | 0.159 (474,672 vox) | 0.157 (473,938 vox) | 0.014 (75,835 vox) |
| 26-connected components / largest share | 116 / **0.999** | 152 / **0.998** | 4,892 / 0.002 |
| largest component bbox (z, y, x) | 1–198, 1–155, 1–198 | 1–198, 1–155, 1–198 | (127 vox speckles) |
| thickness (EDT max / median) | 3.3 / 1.0 | 3.2 / 1.0 | 1.4 / 1.0 |
| raw grey (uint16) mean ± sd — unlabelled fg | 31,485 ± 3,827 | 31,863 ± 4,433 | 32,550 ± 2,477 |
| — labelled grains | 47,368 ± 4,716 | 47,279 ± 4,874 | 46,242 ± 4,483 |
| — non-foreground | 9,062 ± 7,159 | 9,096 ± 7,187 | 22,753 ± 8,095 |

- **Geometry.** A single, thin, continuous sheet (≤ ~6 voxels thick) that spans the full z range and runs diagonally
  across the crop in (y, x). In the z = 100 slice it is an arc from (x ≈ 0, y ≈ 150) to (x ≈ 200, y ≈ 10). In the y = 100
  and x = 100 slices it is a vertical band at x ≈ 105 and y ≈ 105. The 26-connected bounding box stops one voxel short
  of the faces only because of the closing / fill-holes in e13's foreground. In the image the sheet reaches the z = 0,
  z = 199, y = 0 and x = 199 faces: it cuts off the **high-y / high-x corner** of the crop.
- **Intensity.** It sits between grains and air (≈ 31.5k against 47.3k and ≈ 6–9k), as expected for a polymer membrane.
- **Beyond it** is exterior air. GT occupancy falls from 0.68 at y ≈ 0 to 0 for y ≥ 160, and from 0.49 to 0.07 along
  x. The non-foreground is far darker on crop05 (mean 9.1k) than on crop04 (22.8k, pores only), and it pulls the Otsu
  threshold down (119 against 138 on the robust uint8 scale).
- Both scans show the same geometry, as expected: same crop origin, consecutive scans.
- **Renders.** Green = valid GT boundaries, orange = excluded GT, magenta = Otsu foreground without a label.
  - crop05: packed elliptical grains (each with a dark central spot) fill the low-y / low-x part. A magenta membrane
    band separates them from a uniformly black exterior that covers roughly the other half of every slice.
  - crop04: grains everywhere, almost no magenta.

## (b) Stage-by-stage comparison with crop04 (frozen outputs)

| stage | quantity | 00_crop04 | 00_crop05 | 01_crop05 |
|---|---|---|---|---|
| 0 scale / openrec | r_est → openrec radius | 10.2 → 15 | 9.0 → 14 | 8.0 → 12 |
| 1 openrec image | fraction exactly 0 | 0.213 | 0.578 | 0.671 |
| | grain body / pore / exterior (median) | 226 / 8 / – | 227 / 0 / 0 | 196 / 0 / 0 |
| **2 flatten λ=1 (selector)** | **bg_ref = median(blur>0)** | **186.1** | **0.9** | **0.9** |
| | blur>0 voxels below 1 grey level | 0.078 | 0.506 | 0.511 |
| | selector grain body / pore / exterior | 165 / 63 / – | **44 / 0 / 7** | **33 / 5 / 5** |
| 3 ToS (true_dual) | tree nodes | 9,333,317 | 8,877,634 | 8,777,339 |
| | MIN-like / MAX-like nodes | 482,164 / 851,152 | 161,441 / 716,192 | 28,996 / 748,342 |
| | MIN height scale p95 | 32 | 3 | 4 |
| | dark_alt_max (p78 of nonzero) | 178 | 46 | 42 |
| | MIN candidates before → after hierarchy_suppress | 148,487 → 1,800 | 72,619 → 1,800 | 10,635 → 1,800 |
| | selected MIN: h / delta / area / dark score (mean) | 0.95 / 0.71 / 0.50 / 0.89 | 1.00 / 0.07 / 0.82 / 0.95 | 0.73 / 0.12 / 0.62 / 0.96 |
| | selected MIN altitude (median) / area (median) | 17 / 4 | 2 / 24 | 2 / 11 |
| 4 background / boundary rule | dark threshold = min(Otsu − 10, p40, 135) on robust(selector) | 110 | 99 | 91 |
| | dark / external (face-connected) / internal fraction | 0.147 / 0.137 / 0.0096 | 0.709 / 0.703 / 0.0052 | 0.765 / 0.761 / 0.0041 |
| | internal-dark components (total; per valid grain, median) | 19,230; 13 | 17,120; **83** | 15,596; **75** |
| | valid grains containing internal dark | 272 / 272 | 129 / 129 | 128 / 128 |
| | legacy seeds kept / rejected as external | 197 / 407 | **0** / 670 | **0** / 149 |
| 5 e10 spatialisation | MIN regions touching a face | 0 | 0 | 0 |
| | MIN exclusive support (region voxels) median / q75 | 4 / 19 | 24 / 34 | 11 / 20 |
| | MIN points in valid grains (distinct grains) | 348 (268) | **36 (32)** | 155 (85) |
| | MIN points in GT background | 1,269 (70 %) | **1,746 (97 %)** | 1,540 (86 %) |
| 6 e11 point-internal consolidation | points in an internal component → distinct components → markers | 511 → 403 → 398 | 53 → 53 → **53** | 219 → 219 → **191** |
| 7 gate | marker support median / q75 | 43 / 49 (off) | **2 / 3 (fires)** | **2 / 3 (fires)** |
| 8 evaluation | ToS / sel F1@0.50 | 0.985 / 0.985 | 0.012 / 0.260 | 0.066 / 0.807 |

### Reading, stage by stage
1. **Openrec is not the culprit.** Grain bodies keep the same level as crop04 (227 / 196 against 226). The difference
   is that the exterior and the pores are clipped to exactly 0 over 58–67 % of the crop (21 % on crop04). This is the
   precondition, not the break.
2. **Flatten λ = 1 is the first stage that breaks.** It is verified by re-running the frozen function: bit-identical
   selector on all three crops. `bg_ref` is meant to be the background level (186 on crop04, ≈ the grain level), but on
   crop05 the median lands in the float Gaussian tails over the zero plateaux (51 % of "positive" voxels are < 1), giving 0.9.
   Without the offset, `vol − blur + 0.9` is a pure high-pass. After clipping at 0 and the robust p0.5–p99.5 rescale,
   grain bodies end up at 44 / 33 and grain cores, pores and air are merged at 0–7. **The grain-core / pore contrast
   that BA-TOS relies on is gone in the selector image.**
3. **ToS selector (true_dual).** With most of the volume on near-zero plateaux:
   - the MIN-like node count drops 3–17×;
   - the MIN height scale collapses (p95 3–4), so `h_score` saturates at ~1 for everything;
   - `delta_score` falls to 0.07–0.12;
   - `dark_alt_max` (p78 of the nonzero selector voxels) drops to 46 / 42.

   The weighted score then mostly ranks area and darkness. The selected minima sit at altitude ~2, i.e. in the
   pore / exterior plateaux: 86–97 % of the MIN points fall outside every grain.
4. **The boundary rule does not swallow internal candidates.** External dark covers 70–76 % of the crop (14 % on
   crop04). It includes the exterior and the pores, and touches every grain on all three crops, but the internal-dark
   mask still exists inside every valid grain. On crop05, however, the internal mask is **fragments of dim grain bodies**
   (75–83 components per grain), not grain cores. That is a downstream effect of step 2: the dark threshold is applied to
   the compressed selector. The legacy ToS seeds (same compressed image) all fall in the external dark (0 kept), so the
   legacy baseline markers are empty too.
5. **Point-internal consolidation.** Only the few MIN points that land in an internal fragment survive (53 / 219 against 511 on
   crop04), one per fragment → 53 / 191 markers for 129 / 128 grains. Their exclusive supports are the tiny regions of
   those near-plateau minima (median 2) → the gate fires on both crops. **The gate is right to fire, but the collapse is
   decided at step 2.**

### Consistency with earlier findings
- `docs/SCALE_ESTIMATOR.md`: a better openrec radius (16) does not recover the ToS route. The radius does not change
  the zero plateaux, so it does not change `bg_ref`.
- Phase 5a group 2: `flatten_lambda_xy` was the second most sensitive parameter (λ = 0.5: −0.637), and
  `flatten_out_percentile` / `openrec_norm_percentile` mattered most on the min crop. Both are consistent with the
  selector's intensity mapping being fragile.

## Suggested test (not run)
A one-line counterfactual in a copy of the flatten stage, without touching frozen code: compute `bg_ref` over the
specimen only. For example use the median over `blur ≥ 1`, which gives 176.7 / 117.4 on crop05 and 198.8 on crop04, or
exclude the face-connected zero plateau. Re-run from the flatten stage downstream on the two crop05 and on crop04 (which
must stay at 0.985). If crop05 recovers, the root cause is confirmed; if not, the next suspect is the openrec clipping
of the pores to 0.
