# Phase 5a — parameter ranking by effect size

Source: `global_sensitivity.csv`. Effect = max |Δ crop-macro mean F1@0.50| over the perturbations (12 crops, seed-ignore). Influential: > 0.01 on either perturbation, either route (BA-TOS_sel primary).

Identity run (whole pipeline re-run from scratch by the sensitivity driver):

| route | n_crops_done | mean_f1_050 | min_f1_050 | delta_mean_f1_050 | delta_min_f1_050 |
|---|---|---|---|---|---|
| BA-TOS_sel | 12 | 0.9063 | 0.2599 | +0.0000 | +0.0000 |
| BA-TOS_ToS | 12 | 0.8238 | 0.0121 | +0.0000 | +0.0000 |

Method parameters swept and complete: 96 / 96; influential: **25** (sel route only: 25).

## Influential parameters

| group | param | frozen_value | sel_max_abs_delta_mean | tos_max_abs_delta_mean | sel_perturbations |
|---|---|---|---|---|---|
| openrec_flatten | openrec_radius_factor | 1.5 | 0.8827 | 0.8155 | x0.5=0.75: -0.8827, x2=3.0: FAIL, x4/3=2.0: +0.0279 |
| openrec_flatten | flatten_lambda_xy | 1.0 | 0.6373 | 0.5874 | x0.5=0.5: -0.6373, x2=2.0: -0.0154 |
| selector | selector_area_preferred_position | 0.42 | 0.6128 | 0.6535 | x0.5=0.21: -0.5123, x2=0.84: -0.6128 |
| selector | selector_w_delta | 0.12 | 0.5123 | 0.5482 | x0.5=0.06: -0.0070, x2=0.24: -0.5123 |
| selector | selector_w_area | 0.18 | 0.4166 | 0.4524 | x0.5=0.09: -0.4166, x2=0.36: -0.0096 |
| selector | selector_area_preferred_width | 0.55 | 0.4166 | 0.4524 | x0.5=0.275: +0.0017, x2=1.1: -0.4166 |
| selector | selector_area_min | 2 | 0.2887 | 0.3097 | -1=1: -0.2887, +1=3: -0.0408 |
| markers | bg_dark_threshold_max | 135 | 0.2291 | 0.2291 | x0.5=68: -0.2291, x2=255: +0.0000 |
| selector | selector_height_scale_percentile | 95 | 0.2164 | 0.1789 | x0.5=47.5: -0.2164, max=100: -0.0162 |
| selector | selector_K | 3600 | 0.1273 | 0.1196 | x0.5=1800: -0.1273, x2=7200: +0.0321 |
| selector | selector_K_min_fraction | 0.5 | 0.1273 | 0.1196 | x0.5=0.25: -0.1273, x2 (all MIN)=1.0: +0.0020 |
| selector | selector_w_dark | 0.18 | 0.1014 | 0.1139 | x0.5=0.09: -0.0059, x2=0.36: -0.1014 |
| gate_watershed | gate_median_limit | 3.0 | 0.0825 | 0.0000 | x0.5=1.5: -0.0825, x2=6.0: +0.0000 |
| gate_watershed | gate_q75_limit | 4.0 | 0.0825 | 0.0000 | x0.5=2.0: -0.0825, x2=8.0: +0.0000 |
| markers | bg_dark_percentile | 40.0 | 0.0585 | 0.0033 | x0.5=20.0: -0.0585, x2=80.0: FAIL |
| markers | edt_nms_radius | 6.0 | 0.0469 | 0.0000 | x0.5=3.0: -0.0469, x2=12.0: +0.0161 |
| selector | selector_delta_norm | 20.0 | 0.0385 | 0.0095 | x0.5=10.0: +0.0385, x2=40.0: -0.0084 |
| selector | selector_area_max | 2200 | 0.0375 | 0.0317 | x0.5=1100: -0.0375, x2=4400: +0.0147 |
| selector | selector_dark_percentile_max | 78.0 | 0.0370 | 0.0084 | x0.5=39.0: +0.0370, x2=100.0: -0.0006 |
| openrec_flatten | openrec_radius_floor | 8 | 0.0363 | 0.0021 | x0.5=4: +0.0000, x2=16: +0.0363 |
| selector | selector_delta_parent_min | 0 | 0.0312 | 0.0013 | +1=1: +0.0000, +2=2: +0.0312 |
| selector | selector_w_height | 0.22 | 0.0288 | 0.0277 | x0.5=0.11: -0.0288, x2=0.44: -0.0024 |
| openrec_flatten | flatten_out_percentile | 0.5 | 0.0260 | 0.0238 | x0.5=0.25: -0.0260, x2=1.0: +0.0227 |
| openrec_flatten | flatten_sigma_z | 0.0 | 0.0162 | 0.0138 | +1=1.0: -0.0162 |
| markers | e11_nms_radius | 6.0 | 0.0108 | 0.0061 | x0.5=3.0: -0.0108, x2=12.0: +0.0046 |

## All parameters (ranked)

| group | param | kind | frozen_value | sel_max_abs_delta_mean | sel_max_abs_delta_min | tos_max_abs_delta_mean | sel_breaks | complete |
|---|---|---|---|---|---|---|---|---|
| openrec_flatten | openrec_radius_factor | scale | 1.5 | 0.8827 | 0.3160 | 0.8155 | True | True |
| openrec_flatten | flatten_lambda_xy | flag | 1.0 | 0.6373 | 0.1646 | 0.5874 | False | True |
| selector | selector_area_preferred_position | flag | 0.42 | 0.6128 | 0.1398 | 0.6535 | False | True |
| selector | selector_w_delta | flag | 0.12 | 0.5123 | 0.2599 | 0.5482 | False | True |
| selector | selector_w_area | flag | 0.18 | 0.4166 | 0.2599 | 0.4524 | False | True |
| selector | selector_area_preferred_width | flag | 0.55 | 0.4166 | 0.1070 | 0.4524 | False | True |
| selector | selector_area_min | flag | 2 | 0.2887 | 0.2599 | 0.3097 | False | True |
| markers | bg_dark_threshold_max | flag | 135 | 0.2291 | 0.0000 | 0.2291 | False | True |
| gate_watershed | eval_exclude_face_touching | eval | True | 0.2266 | 0.1419 | 0.1964 | False | True |
| selector | selector_height_scale_percentile | literal | 95 | 0.2164 | 0.1771 | 0.1789 | False | True |
| selector | selector_K | flag | 3600 | 0.1273 | 0.4778 | 0.1196 | False | True |
| selector | selector_K_min_fraction | literal | 0.5 | 0.1273 | 0.1914 | 0.1196 | False | True |
| selector | selector_w_dark | flag | 0.18 | 0.1014 | 0.1462 | 0.1139 | False | True |
| gate_watershed | gate_median_limit | literal | 3.0 | 0.0825 | 0.2478 | 0.0000 | False | True |
| gate_watershed | gate_q75_limit | literal | 4.0 | 0.0825 | 0.2478 | 0.0000 | False | True |
| markers | bg_dark_percentile | flag | 40.0 | 0.0585 | 0.2333 | 0.0033 | True | True |
| markers | edt_nms_radius | flag | 6.0 | 0.0469 | 0.1122 | 0.0000 | False | True |
| selector | selector_delta_norm | flag | 20.0 | 0.0385 | 0.4601 | 0.0095 | False | True |
| selector | selector_area_max | flag | 2200 | 0.0375 | 0.1319 | 0.0317 | False | True |
| selector | selector_dark_percentile_max | flag | 78.0 | 0.0370 | 0.5030 | 0.0084 | False | True |
| openrec_flatten | openrec_radius_floor | scale | 8 | 0.0363 | 0.4417 | 0.0021 | False | True |
| selector | selector_delta_parent_min | flag | 0 | 0.0312 | 0.4197 | 0.0013 | False | True |
| selector | selector_w_height | flag | 0.22 | 0.0288 | 0.0587 | 0.0277 | False | True |
| openrec_flatten | flatten_out_percentile | literal | 0.5 | 0.0260 | 0.5243 | 0.0238 | False | True |
| openrec_flatten | flatten_sigma_z | flag | 0.0 | 0.0162 | 0.2330 | 0.0138 | False | True |
| markers | e11_nms_radius | flag | 6.0 | 0.0108 | 0.0185 | 0.0061 | False | True |
| selector | selector_area_factor | flag | 22.0 | 0.0054 | 0.0644 | 0.0000 | False | True |
| selector | selector_min_score | flag | 0.35 | 0.0053 | 0.0000 | 0.0002 | False | True |
| openrec_flatten | openrec_norm_percentile | literal | 1 | 0.0043 | 0.1889 | 0.0075 | False | True |
| markers | bg_dark_connectivity | flag | 1 | 0.0041 | 0.0085 | 0.0012 | True | True |
| gate_watershed | ws_distance_norm_percentile | literal | 1.0 | 0.0037 | 0.0000 | 0.0006 | False | True |
| markers | bg_otsu_offset | flag | -10 | 0.0028 | 0.0000 | 0.0028 | False | True |
| gate_watershed | ws_mask_closing_iterations | literal | 1 | 0.0026 | 0.0113 | 0.0068 | True | True |
| gate_watershed | ws_mask_otsu_offset | literal | 0.0 | 0.0015 | 0.0113 | 0.0035 | False | True |
| gate_watershed | ws_connectivity | literal | 1 | 0.0011 | 0.0113 | 0.0035 | False | True |
| gate_watershed | ws_image_norm_percentile | literal | 1.0 | 0.0003 | 0.0000 | 0.0011 | False | True |
| markers | lm_connectivity | flag | 1 | 0.0003 | 0.0000 | 0.0003 | False | True |
| markers | lm_center_radius | flag | 2 | 0.0003 | 0.0000 | 0.0003 | False | True |
| markers | lm_ring_inner | flag | 4 | 0.0003 | 0.0000 | 0.0003 | False | True |
| markers | lm_center_dark_max | flag | 205.0 | 0.0003 | 0.0000 | 0.0003 | False | True |
| markers | lm_nms_min_dist | flag | 7.0 | 0.0003 | 0.0000 | 0.0003 | False | True |
| markers | lm_max_markers | flag | 420 | 0.0003 | 0.0000 | 0.0003 | False | True |
| selector | selector_w_support | flag | 0.3 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_K | flag | 900 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_min_score | flag | 0.35 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_height_factor | flag | 0.012 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_area_factor | flag | 22.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_area_min | flag | 2 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_area_max | flag | 2200 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_dark_percentile_max | flag | 78.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_delta_parent_min | flag | 0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_support_area_min | flag | 20 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_support_area_max | flag | 45000 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_support_max_steps | flag | 16 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_require_max_support | flag | False | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_support_area_ratio_preferred | flag | 18.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_support_area_ratio_width | flag | 80.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_support_contrast_norm | flag | 80.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_delta_norm | flag | 20.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_area_preferred_position | flag | 0.42 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_area_preferred_width | flag | 0.55 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_w_height | flag | 0.22 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_w_delta | flag | 0.12 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_w_area | flag | 0.18 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_w_dark | flag | 0.18 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_w_support | flag | 0.3 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_component_min | flag | 2 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_component_max | flag | 1500 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_cleanup_open_radius | flag | 0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_cleanup_close_radius | flag | 0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | legacy_support_height_ratio | literal | 0.1 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_seed_area_min | flag | 1 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_seed_area_max | flag | 3000 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_max_external_overlap | flag | 0.1 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_min_internal_overlap | flag | 0.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_internal_search_radius | flag | 6 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_center_radius | flag | 2 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_ring_inner | flag | 4 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_ring_outer | flag | 10 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_center_dark_max | flag | 210.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_min_ring_contrast | flag | -5.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_max_ring_external_coverage | flag | 0.7 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_nms_min_dist | flag | 5.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | bg_seed_connectivity | literal | 1 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | lm_area_min | flag | 3 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | lm_area_max | flag | 900 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | lm_max_dim | flag | 30 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | lm_max_elongation | flag | 8.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | lm_min_fill | flag | 0.005 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | lm_ring_outer | flag | 10 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | lm_min_ring_contrast | flag | -2.0 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | lm_score_dark_weight | literal | 0.02 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | e11_edt_priority_weight | flag | 0.2 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | e11_min_priority_bonus | flag | 0.05 | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | e11_polarity | flag | min | 0.0000 | 0.0000 | 0.0000 | False | True |
| markers | edt_mask_min_component | literal | 500 | 0.0000 | 0.0000 | 0.0000 | False | True |
| gate_watershed | ws_mask_min_component | literal | 500 | 0.0000 | 0.0000 | 0.0000 | False | True |
| gate_watershed | eval_min_gt_voxels | eval | 1000 | 0.0000 | 0.0000 | 0.0000 | False | True |
