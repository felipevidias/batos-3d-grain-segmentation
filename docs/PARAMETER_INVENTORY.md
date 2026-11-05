# BA-TOS parameter inventory (frozen E-series pipeline, `batos-paper-v2-frozen`)

Phase 5a, 2026-11-05. **Question:** how many parameters does BA-TOS really have?

**Method.**
- Every frozen script was read in full. Command lines are those of `src/batos_v2/frozen.py`, i.e. of the recovered `run_e13_e11_all_crops.sh` / `run_batos_pipeline.sh`.
- Listed for each script:
  - every algorithmic **CLI** flag: the value passed in the frozen call, or the script default when the call does not pass it;
  - every **hard-coded** constant that can change the output.
- Script paths are in `SCRIPTS` of `src/batos_v2/frozen.py`.

**Legend.**
- **Origin**
  - `cli`: flag passed in the frozen call.
  - `cli-def`: flag not passed, so the script default applies.
  - `hc`: hard-coded literal.
  - `hc-mod`: module-level constant.
- **Status**
  - `live`: the parameter can change the segmentation.
  - `dead`: proven to have no effect on any output consumed downstream; the reason is given.
  - `via R̂`: acts only through the GT-free scale estimate R̂.
  - `n/a`: categorical or structural choice, not a number.
- **Swept**: the name in `configs/sensitivity_global.yaml`. A dash means not swept, with the reason given.
- **Scales with** is a first guess, to be checked against the sensitivity results (§3):
  - R: grain radius.
  - V: crop volume / object count.
  - I: image intensity (u8 grey levels after normalisation).
  - –: dimensionless or structural.

## 0. Scale estimate — `e13_estimate_scale.py` (linux_patches/BA_TOS_E13_E16…/audit)
GT-free. Its outputs are R̂ (median EDT peak), openrec_radius and min_predicted_size. min_predicted_size is **unused** in the E-series, because Table II is computed with the prediction size filter off.

| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 0.1 | u8 normalisation percentiles | hc | 1 / 99 | via R̂ | – | I |
| 0.2 | Otsu foreground test | hc | `>=` | via R̂ | – | – |
| 0.3 | foreground closing SE | hc | ball(1) | via R̂ | – | R |
| 0.4 | EDT local-max filter size | hc | 5 (5³ cube) | via R̂ | – | R |
| 0.5 | peak floor | hc | max(2.0, P75(EDT)) | via R̂ | – | R |
| 0.6 | min peak count before P95 fallback | hc | 10 | via R̂ | – | V |
| 0.7 | peak cap (top-N) | hc | 1500 | via R̂ | – | V |
| 0.8 | radius statistic | hc | median | via R̂ | – | – |
| 0.9 | R̂ floor | hc | 3.0 | via R̂ | – | R |
| 0.10 | **openrec radius factor** | hc | 1.5 (r = max(8, round(1.5 R̂))) | live | openrec_radius_factor | R |
| 0.11 | **openrec radius floor** | hc | 8 | live (binds only if 1.5 R̂ < 8) | openrec_radius_floor | R |
| 0.12 | min_predicted_size factor / floors | hc | 0.10 · V̂, floor 500, V̂ floor 500 | dead (filter off in E-series) | – | R³ |

## 1. Opening by reconstruction — `batos_preprocess_opening_reconstruction_v2.py` (zsh_recovered)
| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 1.1 | --radius | cli | openrec_radius (15 on crop04; from §0) | live | via 0.10 / 0.11 | R |
| 1.2 | input & output u8 percentiles (same function, twice) | hc | 1 / 99 | live | openrec_norm_percentile | I |
| 1.3 | structuring element | hc | Euclidean ball | n/a | – | R |
| 1.4 | operator | hc | erosion → reconstruction by dilation | n/a | – | – |
| 1.5 | flattening offset | hc | + median(opened) | n/a | – | I |

## 2. Gaussian flattening — `src/batos/batos_preprocess.py` (`--lambda-xy 1`)
| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 2.1 | --lambda-xy | cli | 1.0 (default 15) | live | flatten_lambda_xy | R |
| 2.2 | --sigma-z | cli-def | 0.0 (per-slice 2D) | live | flatten_sigma_z | R |
| 2.3 | output percentiles | hc | 0.5 / 99.5 (on voxels > 0) | live | flatten_out_percentile | I |
| 2.4 | input rescale | hc | 0.5 / 99.5 | dead (uint8 input is passed through) | – | – |
| 2.5 | Gaussian truncate / mode | hc (scipy default) | 4.0 / reflect | n/a | – | – |

## 3. Legacy ToS seed (K = 900) — `batos_tos_seed.py` (pilot_branch_413920f)
- Its seed is consumed by batos_background, giving the internal-dark mask. That mask feeds E10 spatialisation, the E11 markers and the baseline markers.
- `seed_raw.tif` is written **after** `cleanup_seed`, so the cleanup and component flags are live here. In true_dual they are dead.

| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 3.1 | --max-selected-nodes (K_legacy) | cli | 900 | live | legacy_K | V |
| 3.2 | --min-score | cli | 0.35 | live | legacy_min_score | – |
| 3.3 | --height-factor | cli | 0.012 | live | legacy_height_factor | – |
| 3.4 | --area-factor | cli | 22 | live | legacy_area_factor | – |
| 3.5 | --area-min / --area-max | cli | 2 / 2200 | live | legacy_area_min / _max | R³ |
| 3.6 | --delta-parent-min | cli | 0 | live | legacy_delta_parent_min | I |
| 3.7 | --dark-percentile-max | cli | 78 | live | legacy_dark_percentile_max | I (percentile) |
| 3.8 | --support-area-min / --support-area-max | cli | 20 / 45000 | live | legacy_support_area_min / _max | R³ |
| 3.9 | --support-max-steps | cli | 16 | live | legacy_support_max_steps | – |
| 3.10 | --component-min / --component-max | cli | 2 / 1500 | live | legacy_component_min / _max | R³ |
| 3.11 | --cleanup-open-radius / --cleanup-close-radius | cli-def | 0 / 0 | live | legacy_cleanup_open/close_radius | R |
| 3.12 | --require-max-support | cli-def | False | live | legacy_require_max_support | – |
| 3.13 | --support-area-ratio-preferred / -width | cli-def | 18 / 80 | live | legacy_support_area_ratio_preferred / _width | – |
| 3.14 | --support-contrast-norm | cli-def | 80 | live | legacy_support_contrast_norm | I |
| 3.15 | --delta-norm | cli-def | 20 | live | legacy_delta_norm | I |
| 3.16 | --area-preferred-position / -width | cli-def | 0.42 / 0.55 | live | legacy_area_preferred_position / _width | – |
| 3.17 | --w-height / -delta / -area / -dark / -support | cli-def | 0.22 / 0.12 / 0.18 / 0.18 / 0.30 | live | legacy_w_* | – |
| 3.18 | support-height saturation | hc | 0.10 · max MAX height | live | legacy_support_height_ratio | – |
| 3.19 | support sub-weights | hc | 0.35 / 0.30 / 0.25 / 0.10 | live | – (4 coupled weights; see note A) | – |
| 3.20 | --flatten, --flatten-sigma-xy/-z | cli-def | False / 35 / 0 | dead (flatten off) | – | – |
| 3.21 | polarity rule, score sort keys, top-K slice, subtree reconstruction, 6-conn components | hc | – | n/a | – | – |

## 4. Background rejection — `src/batos/batos_background.py`
| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 4.1 | --dark-percentile | cli | 40 | live | bg_dark_percentile | I (percentile) |
| 4.2 | --otsu-offset | cli | −10 | live | bg_otsu_offset | I |
| 4.3 | --dark-threshold-max | cli | 135 | live | bg_dark_threshold_max | I |
| 4.4 | --dark-threshold | cli-def | −1 (automatic) | n/a | – | I |
| 4.5 | --dark-connectivity | cli | 1 | live | bg_dark_connectivity | – |
| 4.6 | --seed-area-min / --seed-area-max | cli | 1 / 3000 | live | bg_seed_area_min / _max | R³ |
| 4.7 | --max-external-overlap | cli | 0.10 | live | bg_max_external_overlap | – |
| 4.8 | --min-internal-overlap | cli | 0.00 | live | bg_min_internal_overlap | – |
| 4.9 | --internal-search-radius | cli | 6 | live | bg_internal_search_radius | R |
| 4.10 | --center-radius / --ring-inner / --ring-outer | cli-def | 2 / 4 / 10 | live | bg_center_radius / bg_ring_inner / bg_ring_outer | R |
| 4.11 | --center-dark-max | cli | 210 | live | bg_center_dark_max | I |
| 4.12 | --min-ring-contrast | cli | −5 | live | bg_min_ring_contrast | I |
| 4.13 | --max-ring-external-coverage | cli | 0.70 | live | bg_max_ring_external_coverage | – |
| 4.14 | --nms-min-dist | cli | 5 | live | bg_nms_min_dist | R |
| 4.15 | seed component connectivity | hc | 6-conn (3, 1) | live | bg_seed_connectivity | – |
| 4.16 | threshold rule | hc | min(Otsu + offset, P_dark, max) clipped 1..255 | n/a | – | I |
| 4.17 | rescale percentiles | hc | 0.5 / 99.5 | dead (u8 input passed through) | – | – |

## 5. Legacy marker builder — `src/batos/batos_markers.py`
Its output (the "baseline markers") is only OR-ed into the watershed mask of both routes. It is never a seed source.

| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 5.1 | --connectivity | cli | 1 | live | lm_connectivity | – |
| 5.2 | --area-min / --area-max | cli | 3 / 900 | live | lm_area_min / _max | R³ |
| 5.3 | --max-dim | cli | 30 | live | lm_max_dim | R |
| 5.4 | --max-elongation | cli | 8 | live | lm_max_elongation | – |
| 5.5 | --min-fill | cli | 0.005 | live | lm_min_fill | – |
| 5.6 | --center-radius / --ring-inner / --ring-outer | cli-def | 2 / 4 / 10 | live | lm_center_radius / lm_ring_inner / lm_ring_outer | R |
| 5.7 | --center-dark-max | cli | 205 | live | lm_center_dark_max | I |
| 5.8 | --min-ring-contrast | cli | −2 | live | lm_min_ring_contrast | I |
| 5.9 | --nms-min-dist | cli | 7 | live | lm_nms_min_dist | R |
| 5.10 | --max-markers | cli | 420 | live | lm_max_markers | V |
| 5.11 | score darkness weight | hc | 0.02 | live | lm_score_dark_weight | – |
| 5.12 | score darkness reference | hc | 205.0 (not tied to --center-dark-max) | live | – (coupled to 5.11) | I |

## 6. ToS candidate selector (true dual, K = 3600) — `batos_tos_seed_true_dual.py` (chat_recovered)
**Score.**
- score = (w_h·h + w_δ·d + w_a·a + w_I·I) / (w_h + w_δ + w_a + w_I).
- h = |height| / P95 per polarity; d = Δparent / δ_norm; a = triangular in log-area; I = darkness, or brightness for MAX.
- Eligibility: area range, area_factor × the polarity's mean area, the intensity percentile, Δparent ≥ min, score ≥ τ_s.

**K.** Applied per polarity: MIN gets ⌊K/2⌋ and MAX gets the rest, each with hierarchy suppression (one node per branch). It is applied again globally after the merge.

Only the candidate CSV is consumed downstream; E10 rebuilds the regions from the CSV.

| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 6.1 | --max-selected-nodes (K) | cli | 3600 | live | selector_K | V |
| 6.2 | --min-score (τ_s) | cli | 0.35 | live | selector_min_score | – |
| 6.3 | --w-height / -delta / -area / -dark | cli | 0.22 / 0.12 / 0.18 / 0.18 | live | selector_w_* | – |
| 6.4 | --area-factor | cli | 22 | live | selector_area_factor | – |
| 6.5 | --area-min / --area-max | cli | 2 / 2200 | live | selector_area_min / _max | R³ |
| 6.6 | --dark-percentile-max | cli | 78 (MAX: mirror P22) | live | selector_dark_percentile_max | I (percentile) |
| 6.7 | --delta-norm | cli | 20 | live | selector_delta_norm | I |
| 6.8 | --delta-parent-min | cli | 0 | live | selector_delta_parent_min | I |
| 6.9 | --area-preferred-position / -width | cli | 0.42 / 0.55 | live | selector_area_preferred_* | – |
| 6.10 | height normaliser percentile | hc | P95 (per polarity) | live | selector_height_scale_percentile | – |
| 6.11 | K split MIN : MAX | hc | ⌊K/2⌋ : rest | live | selector_K_min_fraction | – |
| 6.12 | **--w-support** | cli | 0.30 | **dead**: not in the score (`active_weight` sums 4 weights) | selector_w_support (dead check only) | – |
| 6.13 | **--height-factor** | cli | 0.012 | **dead**: never read | – | – |
| 6.14 | **--support-area-min/-max, --support-max-steps, --support-area-ratio-preferred/-width, --support-contrast-norm, --require-max-support** | cli / cli-def | 20 / 45000 / 16 / 18 / 80 / 80 / False | **dead**: `nearest_max_support` is never called | – | – |
| 6.15 | **--cleanup-open/close-radius, --component-min/-max** | cli | 0 / 0 / 2 / 1500 | **dead**: they act only on the seed TIFF, written after the candidate CSV and unused downstream | – | – |
| 6.16 | --flatten, sigmas | cli-def | False | dead | – | – |
| 6.17 | intensity normalisers, triangular area score, clamp01, sort / tie-break keys, hierarchy rule, Higra ToS defaults (padding mean, immersion) | hc | – | n/a | – | – |

## 7. Spatialisation — `e10_spatialize_candidates_no_gt.py`
No tunable parameter. It rebuilds the exact ToS with Higra defaults and assigns each leaf to its nearest selected ancestor. It reports region_voxels and the centroid-nearest point, with 6-connected internal-dark overlap. **0 parameters.**

## 8. E11 ToS marker builder — `e11_build_markers_ablation.py`
| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 8.1 | --polarity | cli | min | live | e11_polarity (→ dual) | – |
| 8.2 | --nms-radius | cli | 6 | live | e11_nms_radius | R |
| 8.3 | --edt-priority-weight | cli | 0.20 | live | e11_edt_priority_weight | – |
| 8.4 | --min-priority-bonus | cli | 0.05 | live | e11_min_priority_bonus | – |
| 8.5 | one marker per internal component, min-max score normalisation, NMS strict `<` | hc | – | n/a | – | – |

## 9. EDT-Otsu fallback markers — `e11_build_control_markers.py --mode otsu_edt`
These are used only where the gate fires (both crop05s). Their count is matched to the ToS marker count.

| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 9.1 | --nms-radius | cli | 6 | live | edt_nms_radius | R |
| 9.2 | mask min component | hc (func default) | 500 | live | edt_mask_min_component | R³ |
| 9.3 | p1/p99 normalisation, Otsu `>`, 6-conn closing ×1, fill holes, peak_local_max(min_distance 1) | hc | – | live / n/a | – (same constants as §10) | I / – |

## 10. Fixed-mask watershed — `e11_run_fixed_mask_watershed.py` (both routes)
| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 10.1 | --outside-policy | cli | drop | n/a | – | – |
| 10.2 | mask min component | hc (func default) | 500 | live | ws_mask_min_component | R³ |
| 10.3 | mask closing iterations (6-conn) | hc | 1 | live | ws_mask_closing_iterations | R |
| 10.4 | Otsu offset | hc (implicit) | 0 (`image_norm > threshold`) | live | ws_mask_otsu_offset | I |
| 10.5 | image normalisation percentiles | hc | 1 / 99 | live | ws_image_norm_percentile | I |
| 10.6 | distance normalisation percentiles | hc | 1 / 99 (in mask) | live | ws_distance_norm_percentile | R |
| 10.7 | watershed connectivity | hc | 1 | live | ws_connectivity | – |
| 10.8 | watershed line / compactness | hc | False / 0 | n/a | – | – |

## 11. Support-collapse gate — `e23c_build_selective_marker_route.py`
| # | Parameter | Origin | Frozen | Status | Swept | Scales with |
|---|---|---|---|---|---|---|
| 11.1 | MEDIAN_LIMIT | hc-mod | 3.0 voxels | live | gate_median_limit | R³ (support size) |
| 11.2 | Q75_LIMIT | hc-mod | 4.0 voxels | live | gate_q75_limit | R³ |
| 11.3 | statistic population / AND rule / `<=` | hc | post-NMS ToS markers / AND / inclusive | n/a | – | – |

## 12. Evaluator (protocol — not a method parameter) — `batos_v2.evaluate_v2` (≡ e09b, Phase 1d)
| # | Parameter | Frozen | Swept |
|---|---|---|---|
| 12.1 | min GT voxels | 1000 | eval_min_gt_voxels |
| 12.2 | exclude face-touching GT | True | eval_exclude_face_touching |
| 12.3 | ignore rule | dominant GT label of the seed region | – (Phase 1d) |
| 12.4 | IoU thresholds | 0.25 / 0.50 / 0.75 | – (reported) |

## Counts

| | Method parameters |
|---|---|
| Numeric parameters exposed or hard-coded and listed above (§0–§11) | **≈ 125** |
| — of which dead (true_dual support/cleanup/height-factor/w-support, flatten flags, unused input rescales, min_predicted_size) | ≈ 22 |
| — of which act only through R̂ (§0.1–0.9) | 9 |
| — categorical/structural (n/a) | not counted |
| **Swept individually** (configs/sensitivity_global.yaml) | 95 method + 2 evaluator + 1 dead-code check |
| Not swept, live | support sub-weights (4, legacy), score reference 205.0, EDT-Otsu normalisation constants (shared form with §10), R̂ internals (via R̂) |

**Note A.** The four legacy support sub-weights enter one weighted sum inside the legacy score. The legacy score is weighted by `--w-support`, which is swept. Its effect is bounded by legacy_w_support ×0.5/×2.

## Sensitivity results
See §3 below; the section is filled in when each group finishes (`results/sensitivity/global_sensitivity.csv`).
