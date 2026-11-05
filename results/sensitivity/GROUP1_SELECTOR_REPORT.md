# Phase 5a — group 1 (selector weights and K): report

- Sweep: `python -I -m batos_v2.sensitivity --identity --groups selector`, 12 crops, one-factor perturbations of the
  true_dual selector (`configs/sensitivity_global.yaml`, group `selector`), re-run from `true_dual` downstream;
  evaluation `evaluate_v2` with markers (rule dominant, seed-ignore, no size filter == Table II protocol).
- 408 crop-jobs = 12 identity + 33 perturbations × 12; **408/408 complete, 0 failed.**
- Code: launched at 8cceddf + uncommitted Phase 5a code (committed unchanged as 4a4f259); resumed twice on 4a4f259
  (power outage 13:43; Claude Code memory-pressure kill 16:37). Partial job folders were deleted and re-run; completed
  jobs were never re-run. Budget 22 GB / 10 workers until 16:37, then 16 GB / 8 workers (runs are deterministic).
- Raw outputs: `C:\BA_TOS_DATA\runs\sensitivity\<param>_<value>\<crop>\{metrics.csv, run.json, commands.log}`.
  Tables: `global_sensitivity.csv`, `ranking.csv`, `RANKING.md` (this folder).

## Identity

The driver's from-scratch re-run reproduces the frozen results on all 12 crops (per-crop F1@0.50 and F1@0.75 identical,
max |Δ| = 0): sel 0.9063 ± 0.2099 (min 0.2599), ToS 0.8238 ± 0.3670 (median 0.985, min 0.0121), F1@0.75 0.8175 →
0.8992; gate fires only on both crop05 (mv 2, qv 3); closest non-triggered 01_crop03 (17, 43).

## Ranking (effect = max |Δ crop-macro mean F1@0.50|, BA-TOS_sel; threshold 0.01)

**14 of 17 selector parameters are influential.** Not influential: `area_factor` (0.005), `min_score` (0.005),
`w_support` (0.000 — confirms the static dead-code finding: `--w-support` is not in the true_dual score).

| param | frozen | worst perturbation (sel Δ mean) | other direction |
|---|---|---|---|
| area_preferred_position | 0.42 | 0.84: −0.613 | 0.21: −0.512 |
| w_delta | 0.12 | 0.24: −0.512 | 0.06: −0.007 |
| w_area | 0.18 | 0.09: −0.417 | 0.36: −0.010 |
| area_preferred_width | 0.55 | 1.1: −0.417 | 0.275: +0.002 |
| area_min | 2 | 1: −0.289 | 3: −0.041 |
| height_scale_percentile | 95 | 47.5: −0.216 | 100: −0.016 |
| K | 3600 | 1800: −0.127 | 7200: +0.032 |
| K_min_fraction | 0.5 | 0.25: −0.127 | 1.0 (all MIN): +0.002 |
| w_dark | 0.18 | 0.36: −0.101 | 0.09: −0.006 |
| delta_norm | 20 | 10: +0.039 | 40: −0.008 |
| area_max | 2200 | 1100: −0.038 | 4400: +0.015 |
| dark_percentile_max | 78 | 39: +0.037 | 100: −0.001 |
| delta_parent_min | 0 | 2: +0.031 | 1: 0.000 |
| w_height | 0.22 | 0.11: −0.029 | 0.44: −0.002 |

Reading:
- The selector's area prior (`area_preferred_position`, `area_preferred_width`, `w_area`, `area_min`) and the
  `w_delta`/`w_area` balance dominate: single x2 / x0.5 changes collapse the crop-macro mean by 0.3–0.6. The frozen
  operating point sits on a narrow ridge for these.
- K behaves as in Table IV (more candidates help, fewer hurt). `K_min_fraction 0.25` equals `K 1800` exactly (same
  number of MIN candidates, 900), consistent with the MIN-only E11 marker route.
- Several perturbations *improve* the mean (K 7200, delta_norm 10, dark_percentile_max 39, delta_parent_min 2,
  area_max 4400), mostly through the min crop (00_crop05: Δmin up to +0.50). These are GT-informed observations on the
  development crops and must NOT be adopted as new defaults without a held-out test (CLAUDE.md §3 rule 4).

## Gate behaviour under perturbation (input for Phase 5b)

The frozen gate (`median(v_i) ≤ 3 ∧ Q75(v_i) ≤ 4` → EDT-Otsu) changes decision in 5 of 33 perturbations:

| perturbation | crop | gate | mv | qv | preds | sel F1@0.50 | ToS F1@0.50 | identity sel |
|---|---|---|---|---|---|---|---|---|
| w_area 0.18→0.36 | 00_crop05 | fires → **silent** | 18 | 18 | 1 | 0.000 | 0.000 | 0.2599 |
| w_delta 0.12→0.06 | 00_crop05 | fires → **silent** | 14.5 | 16.25 | 2 | 0.000 | 0.000 | 0.2599 |
| delta_norm 20→40 | 00_crop05 | fires → **silent** | 11 | 14.5 | 3 | 0.000 | 0.000 | 0.2599 |
| area_min 2→3 | 00_crop05 | fires → **silent** | 4 | 7 | 21 | 0.000 | 0.000 | 0.2599 |
| area_preferred_position 0.42→0.84 | 01_crop03 | silent → **fires** | 2 | 4 | 182 | 0.652 | 0.542 | 0.9501 |

(00_crop05: 129 valid GT; 01_crop03: 213.) 01_crop05 keeps firing (mv 2, qv 3) in every perturbation.

1. **Count collapse escapes the gate (3 cases).** With 1–3 markers left in the crop, median/Q75 of region voxels are
   computed on 1–3 values and land above 3/4, so the fallback is not used and the crop scores 0. The gate tests support
   *size*, not marker *count*.
2. **The gate limits are coupled to `area_min`.** With `area_min = 3` no candidate can have support < 3, so the
   collapsed crop gives mv 4 / qv 7 with 21 predictions and is missed by a near margin. The absolute 3/4 limits
   are therefore only valid at `area_min = 2`.
3. **One new firing works as intended.** `area_preferred_position 0.84` degrades 01_crop03 into support collapse
   (mv 2, qv 4); the gate fires and the fallback lifts sel above ToS (0.652 vs 0.542), although well below the frozen
   0.950.

Phase 5b candidates: add a count test (markers vs the count expected from the GT-free R̂ and the Otsu mask volume),
and express the support limits relative to `area_min` (or as a continuous support-quality score) instead of fixed
3/4. Freeze any new gate on the development crops before touching Test.
