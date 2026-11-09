# Phase 5e step 4 — frozen BA-TOS v3 on the 52 held-out test crops (frozen / fixflat / v3)

Date 2026-11-08. Branch `journal-extension`. Method under test: **batos-paper-v3-scalefree** (tag
`batos-paper-v3-scalefree`, `configs/batos_paper_v3_scalefree.yaml`, SHA-256 of all executed code recorded at freeze in
db6d6fd). It is compared with **batos-paper-v2-frozen** and **batos-paper-v2-fixflat**, whose Phase-1f test outputs
are reused as they are (deterministic frozen code, not re-run). **No parameter was changed for or after this test.** GT was used only
for evaluation and diagnostic columns (GT inscribed radius, GT counts in the forensics).

## 1. Protocol and run record

* Crops: the same 52 held-out crops as Phase 1f (`data/metadata/crops_test.json`, 8,614 valid GT grains): EFR 12,
  EFS/EMR/EMS/EOR/EOS 8 each. The v3 design used only the 12 dev crops plus the synthetic scale test.
* Runs: `python -I scripts/run_test_v3.py`. Before every session it checks that each frozen file (v3 code, patched copies,
  recovered scripts, the driver and the aggregator) still matches its freeze-time SHA-256. It passed in all sessions
  ("v3 frozen config OK, code commit a853e83"). Evaluation: `batos_v2.evaluate_v2`, IoU 0.25/0.50/0.75, size filter
  none / 500 / own (500 s³); modes `markers_dominant` (seed-ignore, **primary**) and `no_markers_fraction050`.
* Aggregation: `python -I scripts/aggregate_test_v3.py` (frozen, unchanged). Supplementary forensics:
  `python -I scripts/forensics_test_v3_seed_density.py` (new, diagnosis only).
* **Interruptions (four sessions, same settings throughout):**
  1. Started 05:03. The machine shut down after 10 crops. The partial `EFRGP01_01_test11` stopped after flatten and was deleted.
  2. Resumed 12:11. Claude Code stopped the session at 14:36 because the machine was low on memory. 27 crops
     were done; the partial `EMRGP01_01_test08` stopped after scale.
  3. Resumed by the user in a separate terminal at 17:19. Stopped after 29 crops; partial `EMSGP01_00_test02`.
  4. Resumed 20:43, END 23:14.
  Every partial folder was deleted (by hand or by the driver's own clean-up) and that crop was rerun from scratch.
  Completed crops were never re-run. The `test_v3.log` file mixes UTF-8 and UTF-16 lines (the session-3 terminal used
  PowerShell redirection).
* Outcome: 52/52 DONE, **0 crashes, 0 empty results** (fixflat had 5 crashes and frozen 1, all scored as empty).
  624 v3 rows were added to `results/all_runs.csv` (`experiment_id test_v3_v3/<crop>`).
* Dry run: before the run finished, the frozen aggregator was exercised on the first 10 crops in a scratch directory with a
  shim. That run hit a missing `frozen_failed` column: the Phase-1f JSONs for the first 28 crops (EFR, EFS, EMR) predate
  that key. On the full 52-crop set the column exists (it is NaN for those crops, and none of them crashed), so the
  frozen aggregator ran unmodified.

## 2. Results — crop-macro F1@0.50 (markers/dominant = seed-ignore), BA-TOS_sel

v3 ToS ≡ sel on every crop: the gate never fires.

| scope | n | frozen | fixflat | **v3** | v3 median | v3 min | # < 0.90 (frozen / fixflat / v3) |
|---|---|---|---|---|---|---|---|
| EFR_internal | 12 | 0.879 ± 0.223 | 0.956 ± 0.048 | **0.985 ± 0.011** | 0.988 | 0.965 | 4 / 2 / 0 |
| EFS | 8 | 0.882 ± 0.155 | 0.968 ± 0.043 | **0.974 ± 0.044** | 0.993 | 0.876 | 3 / 1 / 1 |
| EMR | 8 | 0.303 ± 0.439 | 0.250 ± 0.416 | **0.980 ± 0.027** | 0.991 | 0.936 | 6 / 6 / 0 |
| EMS | 8 | 0.432 ± 0.465 | 0.479 ± 0.474 | **0.986 ± 0.020** | 0.996 | 0.948 | 6 / 5 / 0 |
| EOR | 8 | 0.293 ± 0.362 | 0.238 ± 0.434 | **0.754 ± 0.240** | 0.766 | 0.431 | 8 / 6 / 4 |
| EOS | 8 | 0.265 ± 0.429 | 0.365 ± 0.499 | **0.754 ± 0.214** | 0.750 | 0.486 | 6 / 5 / 5 |
| external (5 specimens) | 40 | 0.435 ± 0.435 | 0.460 ± 0.473 | **0.890 ± 0.178** | 0.983 | 0.431 | 29 / 23 / 10 |
| all | 52 | 0.538 ± 0.437 | 0.575 ± 0.465 | **0.912 ± 0.161** | 0.985 | 0.431 | 33 / 25 / 10 |

ToS route without the gate (all 52 crops): frozen 0.343, fixflat 0.570, v3 0.912. F1@0.75 (all, sel): 0.492 / 0.542 / 0.868. Size filter 500
(all, sel): 0.555 / 0.587 / 0.914; v3 own filter (500 s³): 0.914. No-marker mode (all, sel, no filter): 0.528 / 0.564 / 0.896.

## 3. Results — pooled objects F1@0.50 (separate table; BA-TOS_sel, no filter)

| scope | GT | frozen P / R / F1 | fixflat P / R / F1 | **v3 P / R / F1** |
|---|---|---|---|---|
| EFR_internal | 2476 | 0.946 / 0.913 / 0.929 | 0.955 / 0.979 / 0.967 | **0.986 / 0.985 / 0.985** |
| EFS | 1793 | 0.821 / 0.951 / 0.881 | 0.972 / 0.959 / 0.965 | **0.981 / 0.959 / 0.970** |
| EMR | 1039 | 0.894 / 0.252 / 0.393 | 0.875 / 0.236 / 0.371 | **0.983 / 0.963 / 0.973** |
| EMS | 1070 | 0.695 / 0.385 / 0.495 | 0.909 / 0.391 / 0.546 | **0.990 / 0.977 / 0.983** |
| EOR | 1120 | 0.679 / 0.236 / 0.350 | 0.902 / 0.197 / 0.324 | **0.848 / 0.645 / 0.733** |
| EOS | 1116 | 0.762 / 0.195 / 0.311 | 0.946 / 0.299 / 0.455 | **0.844 / 0.651 / 0.735** |
| external | 6138 | 0.786 / 0.466 / 0.585 | 0.946 / 0.478 / 0.635 | **0.941 / 0.849 / 0.893** |
| all | 8614 | 0.850 / 0.595 / 0.700 | 0.950 / 0.622 / 0.752 | **0.955 / 0.888 / 0.921** |

With size filter 500 (all): 0.722 / 0.760 / 0.922.

## 4. Gate

| specimen | n | frozen fires | fixflat fires | v3 fires | crashes frozen / fixflat / v3 |
|---|---|---|---|---|---|
| EFR | 12 | 4 | 0 | 0 | 0 / 0 / 0 |
| EFS | 8 | 1 | 0 | 0 | 0 / 0 / 0 |
| EMR | 8 | 6 | 6 | 0 | 0 / 0 / 0 |
| EMS | 8 | 5 | 2 | 0 | 0 / 1 / 0 |
| EOR | 8 | 8 | 4 | 0 | 0 / 1 / 0 |
| EOS | 8 | 6 | 1 | 0 | 1 / 3 / 0 |
| all | 52 | 30 | 13 | **0** | 1 / 5 / 0 |

v3 support statistics stay far from the scaled limits: m_v / q_v are 42–51 on EFR and 77–253 / 166–282 on EO, against limits of
3·max(1,s³) and 4·max(1,s³).

## 5. Scale: r_granulo vs GT median inscribed radius (GT diagnostic only)

| specimen | GT r_in | r_granulo (all crops) | ratio | r_e13 median | ratio |
|---|---|---|---|---|---|
| EFR | 10.63 | 10 | 0.94 | 10.12 | 0.95 |
| EFS | 10.63 | 10 | 0.94 | 9.62 | 0.90 |
| EMR | 14.51 | 14 | 0.96 | 7.55 | 0.52 |
| EMS | 14.70 | 14 | 0.95 | 7.97 | 0.54 |
| EOR | 15.30 | 15 | 0.98 | 7.00 | 0.46 |
| EOS | 15.26 | 15 | 0.98 | 7.07 | 0.46 |

Pearson r(r_granulo, GT r_in) = 0.998 (n = 52). r_e13 gives −0.807 (n = 51; one NaN from the frozen crash on
EOS 00_test03, so the aggregator prints `nan`). Caveat: GT r_in barely varies within a specimen (SD ≤ 0.11), so this
correlation is effectively over three size clusters (10.6 / 14.5–14.7 / 15.3). r_granulo is constant within each specimen,
including on the boundary crops (specimen fraction down to 0.79), and its error is ≤ 6 %. **The scale stage is solved
on the test set.**

## 6. Runtime (s per crop, end to end; peak MB of the largest stage)

| specimen | frozen | fixflat | v3 median | v3 max | v3 peak MB | v3 openrec_lowmem |
|---|---|---|---|---|---|---|
| EFR | 338 | 322 | 349 | 527 | 3598 | 0/12 |
| EFS | 325 | 318 | 488 | 508 | 3614 | 0/8 |
| EMR | 249 | 266 | 517 | 567 | 3648 | 8/8 |
| EMS | 273 | 276 | 387 | 427 | 3648 | 8/8 |
| EOR | 247 | 251 | 389 | 412 | 3636 | 8/8 |
| EOS | 251 | 256 | 390 | 399 | 3339 | 8/8 |
| all | 271 | 280 | **397** | 567 | 3648 | 32/52 |

Totals: frozen 4.0 h, fixflat 4.0 h, v3 6.0 h. Every crop with openrec radius 22 (EM/EO) used `openrec_lowmem`. The v3
timings span four sessions on a machine under varying load (one session ended on a low-memory stop), so they are not
contention-controlled. The peak of 3.65 GB is below the 8 GB budget.

## 7. Flagged crops (v3 F1@0.50 < 0.90): 10 of 52

### 7a. Frozen aggregator rule (as pre-registered)

The frozen rule labels all 10 "watershed: under-segmentation / FN" (recall < precision, after its scale / stage-03 /
markers-vs-valid-GT / gate checks pass). Full table: `tables.md` §Flagged and `forensics_stages.csv`. That rule compares
ToS markers with **valid** GT only (e.g. 156 vs 170), and that comparison misses the actual first anomaly.

### 7b. Supplementary forensics: the first anomalous stage is the candidate capacity K (`seed_density.md`)

ToS markers seed every grain, including the face-touching grains that are excluded from scoring. Per grain (valid +
excluded GT, counts used for diagnosis only):

| specimen / class | s | K_MIN cap 1800/s³ | saturated | GT all | candidates / grain | ToS markers / grain | v3 F1 median | flagged |
|---|---|---|---|---|---|---|---|---|
| EFR interior | 1.0 | 1800 | 6/6 | 658 | 2.74 | 0.60 | 0.985 | 0 |
| EFS interior | 1.0 | 1800 | 4/4 | 654 | 2.75 | 0.58 | 0.964 | 1 |
| EMR / EMS interior | 1.4 | 656 | 8/8 | 465 | 1.41 | 0.54–0.56 | 0.960–0.976 | 0 |
| EOR / EOS boundary | 1.5 | 533 | 8/8 | 270–280 | 1.91–1.97 | 0.58 | 0.978–0.983 | 1 |
| **EOR / EOS interior** | 1.5 | 533 | 8/8 | 448–449 | **1.19** | **0.35–0.37** | **0.548–0.556** | **8** |

* **9 of 10 flags have one first stage: candidate capacity.** `true_dual --max-selected-nodes` is scaled v2/s³
  (class N in `src/batos_v3/constants.py`), and the MIN candidate count reaches its cap on all 52 crops. The rule
  assumes the number of grains per 200³ crop falls as s⁻³. On the test set it does not: relative to EFR, interior grain
  counts fall to 0.69–0.70 against s⁻³ = 0.30–0.36. That is because face-touching grains, which are about 60 % of all
  grains, scale roughly as a surface (s⁻²), and partial grains count as well.
  At s = 1.5 the cap leaves 1.19 candidates per grain on interior crops. After boundary rejection and NMS that is
  0.35 markers per grain instead of about 0.58. About a third of the grains get no seed, so their neighbours absorb them
  (merges). Precision is 0.63–0.84 and recall 0.33–0.63, with ~105 evaluable predictions per 166 valid GT. EM
  (1.41 candidates per grain) is just above the point where this breaks. EO boundary crops (fewer grains, 1.9 per
  grain) are not affected.
  Flagged: EORGP01_00_test01 0.558, EORGP01_00_test03 0.596, EORGP01_01_test05 0.431, EORGP01_01_test07 0.553,
  EOSGP01_00_test01 0.724, EOSGP01_00_test03 0.599, EOSGP01_01_test05 0.486, EOSGP01_01_test07 0.498, and
  EOSGP01_01_test08 0.776 (boundary, borderline: 1.66 candidates per grain and 0.53 markers per grain; this is the one
  crop where v3 is **worse than fixflat**, 0.982).
* **EFSGP01_01_test05 (0.876)** has no capacity problem (2.74 candidates per grain). v3, fixflat and frozen all score
  0.874–0.876 there, so it is a crop-intrinsic watershed under-segmentation (P 0.92 / R 0.83) and not a v3 regression.
* The scale (r_granulo 15 vs GT 15.3), stage-03 (≥ 0.94 of seeds kept) and the gate are all normal on every flagged crop.

**Not acted on.** This diagnosis comes from the test set. Changing how K scales (for example s⁻² or a measured density) is
a method change. It would need dev-only design, a re-freeze and a new held-out set (CLAUDE.md §3.4; §4 Phase 5).

## 8. Summary

* v3 removes every crash and every gate firing. It solves the scale estimate (r_granulo within 6 % of GT r_in on all six
  specimens, where e13 is off by about 2× on EM/EO). On the 52 held-out crops, crop-macro F1@0.50 goes from 0.538 (frozen)
  and 0.575 (fixflat) to **0.912**; pooled F1 goes from 0.700 / 0.752 to **0.921**. EFR, EFS, EMR and EMS are all
  ≥ 0.974 crop-macro.
* EO interior crops remain the failure mode (crop-macro EOR 0.754, EOS 0.754). The candidate cap scaled by s⁻³ is the first
  stage to fail. One crop (EOSGP01_01_test08) regresses against fixflat.

Files: `tables.md` (all scopes, filters and modes), `crop_macro.csv`, `pooled.csv`, `per_crop_primary.csv`,
`per_crop_all_modes.csv`, `crop_info.csv`, `scale.csv`, `runtime.csv`, `forensics_stages.csv`, `flags.csv`,
`seed_density.{csv,md}`, `seed_density_by_class.csv`.
