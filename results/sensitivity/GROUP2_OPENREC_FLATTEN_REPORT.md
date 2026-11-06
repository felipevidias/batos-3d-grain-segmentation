# Phase 5a — group 2 (openrec / flatten): report

- Sweep: `python -I -m batos_v2.sensitivity --groups openrec_flatten`, 12 crops, one-factor perturbations re-run from
  the first affected stage (openrec or flatten) downstream; evaluation `evaluate_v2` with markers (rule dominant,
  seed-ignore, no size filter == Table II protocol). Deltas are against the frozen Table II values (identity verified
  12/12 in group 1).
- Runs: main sweep 2026-11-05 19:34–22:28 (16 GB / 8 workers, commit efc2a3c); radius_factor 2.0 added in cd947a8 and run
  with `--values 2.0` (2 workers, 9/12 before a Claude Code memory-pressure kill; the remaining 3 with 1 worker in a
  separate terminal, 2026-11-06 00:30–00:53). 6 parameters, 12 perturbations, 144 crop-jobs: **133 complete, 10
  infeasible, 1 not run** (all in radius_factor 3.0, below).
- Tables: `global_sensitivity.csv`, `ranking.csv`, `RANKING.md` (this folder; at the time of this commit they also hold
  partial group-3 rows, marked `complete = False`).

## radius_factor 3.0: infeasible (resource limit, not a segmentation failure)

At factor 3.0 the openrec radius becomes 24–31 voxels (frozen 12–16). The frozen openrec script
(`batos_preprocess_opening_reconstruction_v2.py`, `opening_by_reconstruction`, line 35: scipy grey erosion with a ball
footprint) raises `MemoryError` at r ≥ 27 on the 32 GB machine. This happens deterministically, within ~1 s, even with a
single worker. Of the r = 24 crops, 01_crop03 completed (openrec alone 732 s) and 01_crop05 was not run, because it was the
job running when the memory guard killed the shell. Making it feasible would require changing frozen code (rule 1).
The ranking script prints this perturbation as `FAIL`/`breaks`; read it as **infeasible**. The effect size of
radius_factor below therefore comes from 0.75 and 2.0 only.

## Ranking (effect = max |Δ crop-macro mean F1@0.50|, BA-TOS_sel; threshold 0.01)

**5 of 6 parameters are influential.** Not influential: `openrec_norm_percentile` (0.004).

| param | frozen | perturbations (sel Δ mean F1@0.50) | ToS max |Δ| | gate fires on (sel) |
|---|---|---|---|---|
| openrec_radius_factor | 1.5 | 0.75: **−0.883**; 2.0: +0.028; 3.0: infeasible | 0.816 | 10 / 2 / – |
| flatten_lambda_xy | 1.0 | 0.5: **−0.637**; 2.0: −0.015 | 0.587 | 2 / 1 |
| openrec_radius_floor | 8 | 4: 0.000 (non-binding); 16: +0.036 | 0.002 | 2 / 2 |
| flatten_out_percentile | 0.5 | 0.25: −0.026; 1.0: +0.023 | 0.024 | 2 / 2 |
| flatten_sigma_z | 0.0 | 1.0: −0.016 | 0.014 | 1 |
| openrec_norm_percentile | 1 | 0.5: −0.004; 2.0: +0.004 | 0.008 | 2 / 2 |

Reading:
- **The openrec radius is the most sensitive single setting so far.** At factor 0.75 the radius 0.75·R̂ falls below
  the floor, so every crop runs at r = 8 (frozen 12–16). Background subtraction then removes grain-scale structure,
  the marker stage collapses (1–92 predictions against 128–272 valid GT), and the crop-macro mean falls to 0.024 (sel)
  / 0.008 (ToS). Going up is mildly positive (2.0: +0.028) until it becomes infeasible (3.0).
- **Halving the flattening λ (0.5) is almost as harmful (−0.64)**; doubling it is neutral for sel and +0.040 for ToS.
- `openrec_radius_floor` 4 is identical to frozen by construction: round(1.5·R̂) ≥ 12 > 4 on every crop, so the floor
  never binds. At 16 it raises the radius on the crops below 16 and gives +0.036 (min crop 00_crop05 0.26 → 0.70).
  This is a GT-informed observation on development crops; do NOT adopt it as a default without held-out test data
  (CLAUDE.md §3 rule 4).

## Gate behaviour under perturbation (input for Phase 5b)

| perturbation | crop(s) | gate | mv / qv | preds | sel F1@0.50 (identity) |
|---|---|---|---|---|---|
| radius_factor 0.75 | 00_crop02, 00_crop04, 00_crop06, 00_crop07, 01_crop02, 01_crop04, 01_crop06, 01_crop07 | silent → **fires** (8 crops) | 2–3 / 2–3.5 | 7–50 | 0.000–0.068 (0.967–0.991) |
| lambda_xy 2.0 | 00_crop05 | fires → **silent** | 4 / 7 | 897 | 0.204 (0.260) |
| sigma_z 1.0 | 00_crop05 | fires → **silent** | 3 / 5 | 337 | 0.027 (0.260) |

1. **The gate detects the collapse, but the fallback cannot repair it.** At radius_factor 0.75 the gate fires on 10
   of 12 crops (8 new plus both crop05), which is a correct detection. The EDT-Otsu fallback is built on the same
   degraded openrec image, though, so sel stays at 0.00–0.07. The gate protects against selector collapse, not
   against preprocessing collapse.
2. **More near-misses on 00_crop05.** With λ 2.0 (mv 4, qv 7) and σz 1.0 (mv 3, qv 5 — misses only on Q75 by 1 voxel)
   the gate goes silent. In group 1 this happened with area_min 3 (4 / 7). The fixed 3 / 4 limits leave no margin on
   this crop. This supports a continuous support-quality score, or limits relative to `area_min`, in Phase 5b.
3. 00_crop03 and 01_crop03 at radius_factor 0.75 (92 and 4 predictions, F1 0) do not trigger the gate. The second is
   another marker-count collapse that escapes the gate, as in group 1.
