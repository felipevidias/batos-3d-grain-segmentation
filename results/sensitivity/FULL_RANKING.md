# Phase 5a — full global sensitivity ranking (groups 1–4, frozen `batos-paper-v2`, 12 dev crops)

Source: `global_sensitivity.csv` (all four groups), ranked by `scripts/rank_sensitivity.py` → `ranking.csv`, `RANKING.md`
(regenerated 2026-11-09 from the complete CSV). One parameter perturbed at a time on all 12 EFR dev crops, evaluate_v2
with markers (seed-ignore), size filter off. Effect = max over the perturbations of |Δ crop-macro mean F1@0.50|;
influential = effect > 0.01 on either route (BA-TOS_sel primary). Sensitivity runs are not in `results/all_runs.csv`
(CLAUDE.md rule 6 exception); per-job provenance is in `C:\BA_TOS_DATA\runs\sensitivity\<param>_<value>\<crop>\`.

Identity (driver re-runs the whole frozen pipeline): sel 0.9063 (min 0.2599), ToS 0.8238 (min 0.0121), Δ = 0 on both.

## Completeness

| group | parameters | perturbations | crop-jobs done | failed crop-jobs |
|---|---|---|---|---|
| 1 selector | 17 | 33 | 396 / 396 | 0 |
| 2 openrec_flatten | 6 | 12 | 133 / 144 | 11 (radius_factor 3.0: 10 MemoryError at r ≥ 27, 1 not run) |
| 3 markers | 65 | 125 | 1497 / 1500 | 3 (see below) |
| 4 gate_watershed | 10 (8 method + 2 eval) | 19 | 216 / 228 | 12 (see below) |

All configured perturbations of groups 3 and 4 were attempted. **96 / 96 method parameters have a complete ranking
row; 25 are influential** (+1 evaluation-protocol parameter, `eval_exclude_face_touching`, reported separately).

### Failed perturbations — classification

| perturbation | crops failed | error | classification |
|---|---|---|---|
| `openrec_radius_factor` 3.0 | 10 + 1 not run | `MemoryError` in frozen openrec (scipy grey erosion, ball r ≥ 27) | infeasible on this machine (32 GB); not a method break |
| `bg_dark_percentile` 80 | both crop05 | e11: "marker volume empty after applying the mask" | **real pipeline break** |
| `bg_dark_connectivity` 3 | EFRGP01_01_crop05 | e11: same error (154 ToS markers, all dropped) | **real pipeline break** (crop05 only) |
| `ws_mask_closing_iterations` 0 ("no closing") | all 12 | e11: "Otsu bright mask empty" | **sweep artifact, invalid point.** scipy `binary_closing(iterations=0)` means *repeat until convergence*, not *no closing*: dilation fills the volume and erosion (border 0) then empties it. "No closing" needs the call removed, not `iterations=0`. The `+1 = 2` point is valid (Δ −0.0026). Not re-run. |

## Influential parameters (effect > 0.01), ranked by BA-TOS_sel

| group | parameter | frozen | sel effect | ToS effect | perturbations (sel Δ) |
|---|---|---|---|---|---|
| openrec_flatten | openrec_radius_factor | 1.5 | 0.8827 | 0.8155 | 0.75: −0.8827, 3.0: FAIL (memory), 2.0: +0.0279 |
| openrec_flatten | flatten_lambda_xy | 1.0 | 0.6373 | 0.5874 | 0.5: −0.6373, 2.0: −0.0154 |
| selector | area_preferred_position | 0.42 | 0.6128 | 0.6535 | 0.21: −0.5123, 0.84: −0.6128 |
| selector | w_delta | 0.12 | 0.5123 | 0.5482 | 0.06: −0.0070, 0.24: −0.5123 |
| selector | w_area | 0.18 | 0.4166 | 0.4524 | 0.09: −0.4166, 0.36: −0.0096 |
| selector | area_preferred_width | 0.55 | 0.4166 | 0.4524 | 0.275: +0.0017, 1.1: −0.4166 |
| selector | area_min | 2 | 0.2887 | 0.3097 | 1: −0.2887, 3: −0.0408 |
| markers | bg_dark_threshold_max | 135 | 0.2291 | 0.2291 | 68: −0.2291, 255: +0.0000 |
| selector | height_scale_percentile | 95 | 0.2164 | 0.1789 | 47.5: −0.2164, 100: −0.0162 |
| selector | K | 3600 | 0.1273 | 0.1196 | 1800: −0.1273, 7200: +0.0321 |
| selector | K_min_fraction | 0.5 | 0.1273 | 0.1196 | 0.25: −0.1273, 1.0: +0.0020 |
| selector | w_dark | 0.18 | 0.1014 | 0.1139 | 0.09: −0.0059, 0.36: −0.1014 |
| gate_watershed | gate_median_limit | 3.0 | 0.0825 | 0.0000 | 1.5: −0.0825, 6.0: +0.0000 |
| gate_watershed | gate_q75_limit | 4.0 | 0.0825 | 0.0000 | 2.0: −0.0825, 8.0: +0.0000 |
| markers | bg_dark_percentile | 40 | 0.0585 | 0.0033 | 20: −0.0585, 80: FAIL (crop05) |
| markers | edt_nms_radius | 6 | 0.0469 | 0.0000 | 3: −0.0469, 12: +0.0161 |
| selector | delta_norm | 20 | 0.0385 | 0.0095 | 10: +0.0385, 40: −0.0084 |
| selector | area_max | 2200 | 0.0375 | 0.0317 | 1100: −0.0375, 4400: +0.0147 |
| selector | dark_percentile_max | 78 | 0.0370 | 0.0084 | 39: +0.0370, 100: −0.0006 |
| openrec_flatten | openrec_radius_floor | 8 | 0.0363 | 0.0021 | 4: +0.0000, 16: +0.0363 |
| selector | delta_parent_min | 0 | 0.0312 | 0.0013 | 1: +0.0000, 2: +0.0312 |
| selector | w_height | 0.22 | 0.0288 | 0.0277 | 0.11: −0.0288, 0.44: −0.0024 |
| openrec_flatten | flatten_out_percentile | 0.5 | 0.0260 | 0.0238 | 0.25: −0.0260, 1.0: +0.0227 |
| openrec_flatten | flatten_sigma_z | 0 | 0.0162 | 0.0138 | 1: −0.0162 |
| markers | e11_nms_radius | 6 | 0.0108 | 0.0061 | 3: −0.0108, 12: +0.0046 |
| *eval protocol* | eval_exclude_face_touching | true | 0.2266 | 0.1964 | false: −0.2266 (protocol, not a method parameter) |

Per group: selector 14/17 influential, openrec_flatten 5/6, markers 4/65, gate_watershed 3/8 method parameters.

## Groups 3–4 findings

- **The marker stage is almost insensitive.** 61 of 65 marker parameters are below 0.01. All 25 `legacy_*` (tos_seed)
  parameters have Δ = 0.0000 on all 12 crops on both routes. All `lm_*` (markers_legacy) parameters are ≤ 0.0003.
  Whether the tos_seed output reaches the scored routes at all was not checked (no static trace was done), so
  "no measurable effect" is the claim, not "dead code".
- **The influential marker parameter is the stage-03 dark-threshold cap** `bg_dark_threshold_max` (135 → 68: −0.229).
  This is the same `min(Otsu−10, p40, 135)` rule whose p40 branch caused the EM/EO collapse on test set 1
  (`results/test_protocol/REPORT.md`). `bg_dark_percentile` (the p40 term) is the only marker parameter that breaks
  the pipeline (80 → empty markers on both crop05).
- **The watershed and mask stage is robust.** Mask Otsu offset ±0.04, normalisation percentiles, connectivity
  2/3 and extra closing all give |Δ| ≤ 0.004.
- **Gate limits.** Halving either limit stops the gate on both crop05 (gate fires 0/12, sel → 0.8238 = ToS, −0.0825).
  Doubling them changes nothing (still 2/12). The 3/4 limits are not on a knife edge from above. From below they sit
  just above the crop05 values (mv 2, qv 3).
- **The evaluation protocol matters more than most method parameters.** Including face-touching GT gives
  −0.227. This supports reporting the exclusion rule explicitly in the journal.

## Not done / open

- `ws_mask_closing_iterations` "no closing": invalid sweep point (see above). It needs a patch that removes the
  closing call; 12 ws-stage crop-jobs, cheap. Not run without the user's decision.
- `openrec_radius_factor` 3.0: infeasible at 32 GB with the frozen openrec. The low-memory openrec
  (`batos_v2/openrec_lowmem.py`) could run it, but that would change the implementation under test.
