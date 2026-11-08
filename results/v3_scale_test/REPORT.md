# Phase 5e step 3 — synthetic scale test: fixflat vs BA-TOS v3 (dev + synthetic only)

**Data.** 12 EFR dev crops at f = 1.0 and their rescaled copies at f = 0.7 / 1.4 / 2.0 (36 synthetic crops,
`scripts/make_synthetic_scale.py`, `data/metadata/crops_synth_scale.json`): 48 crop evaluations, all complete. No test
crop was read. Valid GT ≥ round(1000 f³) voxels, so the same physical grains are scored at every factor.

**Methods.** fixflat = `batos-paper-v2-fixflat` (e13 scale, v2 constants) via `batos_v3.pipeline --preset fixflat`
(bit-identical to the tagged method; f = 1.0 reuses the tagged runs). v3 = `batos_v3.pipeline --preset v3`
(r_granulo scale, constants in grain units, `docs/V3_SCALE_FREE.md`). No parameter changed during step 3.

**Commands.** `python -I scripts/run_v3_scale_test.py` (paused 2026-11-07 14:19 at 25/48 for a PC shutdown; the one
partial run dir, `fixflat/x2/EFRGP01_00_crop03_x2`, was deleted and the remaining 23 crops were rerun on 2026-11-07/08,
23:53 → 05:01), then `python -I scripts/aggregate_v3_scale_test.py`. Tables: `tables.md`; CSVs in this folder.
Evaluation: evaluate_v2, markers/dominant (seed-ignore) primary; the ToS and sel routes are identical on every crop
(the gate never fires).

## Crop-macro F1@0.50 (BA-TOS_sel, sample SD)

| f | fixflat, no filter | v3, no filter | fixflat, filter 500 | v3, filter 500 s³ |
|---|---|---|---|---|
| 0.7 | 0.999 ± 0.001 (min 0.996) | 0.975 ± 0.014 (min 0.946) | 1.000 ± 0.001 | 0.995 ± 0.004 (min 0.987) |
| 1.0 | 0.976 ± 0.023 (min 0.913) | 0.989 ± 0.008 (min 0.974) | 0.983 ± 0.014 | 0.990 ± 0.008 |
| 1.4 | 0.750 ± 0.372 (min 0.000) | 0.995 ± 0.007 (min 0.982) | 0.824 ± 0.385 | 1.000 ± 0.000 |
| 2.0 | 0.239 ± 0.369 (min 0.000) | 0.938 ± 0.138 (min 0.621) | 0.383 ± 0.471 | 1.000 ± 0.000 |

## Pooled objects F1@0.50 (separate table)

| f | GT | fixflat no filter (P / R) | v3 no filter (P / R) | fixflat 500 | v3 500 s³ |
|---|---|---|---|---|---|
| 0.7 | 8939 | 0.999 (0.999 / 1.000) | 0.978 (0.962 / 0.994) | 1.000 | 0.994 |
| 1.0 | 2808 | 0.979 (0.979 / 0.978) | 0.988 (0.991 / 0.985) | 0.984 | 0.989 |
| 1.4 | 877 | 0.819 (0.834 / 0.805) | 0.995 (0.990 / 1.000) | 0.889 | 1.000 |
| 2.0 | 188 | 0.277 (0.234 / 0.340) | 0.947 (0.900 / 1.000) | 0.487 | 1.000 |

## Scale, gate, failures, runtime

| f | GT r_in (median) | r_granulo (all 12 crops) | e13 r_est (median) | fixflat crashes | v3 empty / crash | gate fires (fixflat / v3) | v3 s/crop median (max) | v3 peak MB |
|---|---|---|---|---|---|---|---|---|
| 0.7 | 7.35 | 7 | 7.78 | 0 | 0 / 0 | 0 / 0 | 351 (369) | 3378 |
| 1.0 | 10.63 | 10 | 10.17 | 0 | 0 / 0 | 0 / 0 | 397 (485) | 3355 |
| 1.4 | 14.44 | 14 | 10.00 | 2 | 0 / 0 | 0 / 0 | 444 (550) | 3608 |
| 2.0 | 20.52 | 20–21 | 12.60 | 6 | 0 / 0 | 0 / 0 | 644 (730) | 3610 |

Runtimes are contention-inflated: the second half ran alongside the user's sensitivity sweep (6 workers).
v3 openrec uses `openrec_lowmem` at f ≥ 1.4 (radius 21 / 30–32).

## Reading

* **Scale.** r_granulo follows the true scale exactly (ratio to GT r_in 0.95–0.98 at every factor). e13 saturates
  (10.0 at f 1.4, 12.6 at f 2.0), which is what breaks fixflat at large scales: 8 "no MIN candidate inside internal_dark"
  crashes and collapsed F1 on most of the remaining crops.
* **v3 at f ≥ 1.0:** recall 1.000 at f 1.4 and 2.0. Every v3 FP is a fragment below the scaled size filter.
  With filter 500 s³, F1 is 1.000 on all 24 large-scale crops. The two cells below 0.90 (both crop05 at f 2.0:
  0.667 and 0.621, no filter) are 9-GT-grain crops with 9 and 11 sub-filter fragments.
* **v3 at f 0.7, no filter:** 0.975 against fixflat's 0.999. All of the gap is small FP fragments (351 pooled vs 13;
  recall 0.994). With the method's own filter (171 voxels): 0.995 vs 1.000. This is the one place v3 is behind fixflat.
  It is reported, not tuned.
* **Step-3 criterion** (v3 matches v2 at f = 1 and does not collapse at other scales): met. No design iteration was
  made. v3 is frozen as is for the test protocol (step 4).
