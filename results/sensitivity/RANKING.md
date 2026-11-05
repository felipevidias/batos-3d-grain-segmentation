# Phase 5a — parameter ranking by effect size

Source: `global_sensitivity.csv`. Effect = max |Δ crop-macro mean F1@0.50| over the perturbations (12 crops, seed-ignore). Influential: > 0.01 on either perturbation, either route (BA-TOS_sel primary).

Identity run (whole pipeline re-run from scratch by the sensitivity driver):

| route | n_crops_done | mean_f1_050 | min_f1_050 | delta_mean_f1_050 | delta_min_f1_050 |
|---|---|---|---|---|---|
| BA-TOS_sel | 12 | 0.9063 | 0.2599 | +0.0000 | +0.0000 |
| BA-TOS_ToS | 12 | 0.8238 | 0.0121 | +0.0000 | +0.0000 |

Method parameters swept and complete: 17 / 18; influential: **14** (sel route only: 14).

## Influential parameters

| group | param | frozen_value | sel_max_abs_delta_mean | tos_max_abs_delta_mean | sel_perturbations |
|---|---|---|---|---|---|
| selector | selector_area_preferred_position | 0.4200 | 0.6128 | 0.6535 | x0.5=0.21: -0.5123, x2=0.84: -0.6128 |
| selector | selector_w_delta | 0.1200 | 0.5123 | 0.5482 | x0.5=0.06: -0.0070, x2=0.24: -0.5123 |
| selector | selector_w_area | 0.1800 | 0.4166 | 0.4524 | x0.5=0.09: -0.4166, x2=0.36: -0.0096 |
| selector | selector_area_preferred_width | 0.5500 | 0.4166 | 0.4524 | x0.5=0.275: +0.0017, x2=1.1: -0.4166 |
| selector | selector_area_min | 2.0000 | 0.2887 | 0.3097 | -1=1: -0.2887, +1=3: -0.0408 |
| selector | selector_height_scale_percentile | 95.0000 | 0.2164 | 0.1789 | x0.5=47.5: -0.2164, max=100: -0.0162 |
| selector | selector_K | 3600.0000 | 0.1273 | 0.1196 | x0.5=1800: -0.1273, x2=7200: +0.0321 |
| selector | selector_K_min_fraction | 0.5000 | 0.1273 | 0.1196 | x0.5=0.25: -0.1273, x2 (all MIN)=1.0: +0.0020 |
| selector | selector_w_dark | 0.1800 | 0.1014 | 0.1139 | x0.5=0.09: -0.0059, x2=0.36: -0.1014 |
| selector | selector_delta_norm | 20.0000 | 0.0385 | 0.0095 | x0.5=10.0: +0.0385, x2=40.0: -0.0084 |
| selector | selector_area_max | 2200.0000 | 0.0375 | 0.0317 | x0.5=1100: -0.0375, x2=4400: +0.0147 |
| selector | selector_dark_percentile_max | 78.0000 | 0.0370 | 0.0084 | x0.5=39.0: +0.0370, x2=100.0: -0.0006 |
| selector | selector_delta_parent_min | 0.0000 | 0.0312 | 0.0013 | +1=1: +0.0000, +2=2: +0.0312 |
| selector | selector_w_height | 0.2200 | 0.0288 | 0.0277 | x0.5=0.11: -0.0288, x2=0.44: -0.0024 |

## All parameters (ranked)

| group | param | kind | frozen_value | sel_max_abs_delta_mean | sel_max_abs_delta_min | tos_max_abs_delta_mean | sel_breaks | complete |
|---|---|---|---|---|---|---|---|---|
| selector | selector_area_preferred_position | flag | 0.4200 | 0.6128 | 0.1398 | 0.6535 | False | True |
| selector | selector_w_delta | flag | 0.1200 | 0.5123 | 0.2599 | 0.5482 | False | True |
| selector | selector_w_area | flag | 0.1800 | 0.4166 | 0.2599 | 0.4524 | False | True |
| selector | selector_area_preferred_width | flag | 0.5500 | 0.4166 | 0.1070 | 0.4524 | False | True |
| selector | selector_area_min | flag | 2.0000 | 0.2887 | 0.2599 | 0.3097 | False | True |
| selector | selector_height_scale_percentile | literal | 95.0000 | 0.2164 | 0.1771 | 0.1789 | False | True |
| selector | selector_K | flag | 3600.0000 | 0.1273 | 0.4778 | 0.1196 | False | True |
| selector | selector_K_min_fraction | literal | 0.5000 | 0.1273 | 0.1914 | 0.1196 | False | True |
| selector | selector_w_dark | flag | 0.1800 | 0.1014 | 0.1462 | 0.1139 | False | True |
| selector | selector_delta_norm | flag | 20.0000 | 0.0385 | 0.4601 | 0.0095 | False | True |
| selector | selector_area_max | flag | 2200.0000 | 0.0375 | 0.1319 | 0.0317 | False | True |
| selector | selector_dark_percentile_max | flag | 78.0000 | 0.0370 | 0.5030 | 0.0084 | False | True |
| selector | selector_delta_parent_min | flag | 0.0000 | 0.0312 | 0.4197 | 0.0013 | False | True |
| selector | selector_w_height | flag | 0.2200 | 0.0288 | 0.0587 | 0.0277 | False | True |
| selector | selector_area_factor | flag | 22.0000 | 0.0054 | 0.0644 | 0.0000 | False | True |
| selector | selector_min_score | flag | 0.3500 | 0.0053 | 0.0000 | 0.0002 | False | True |
| selector | selector_w_support | flag | 0.3000 | 0.0000 | 0.0000 | 0.0000 | False | True |
| gate_watershed | gate_median_limit | literal | 3.0000 | nan | nan | nan | False | False |
