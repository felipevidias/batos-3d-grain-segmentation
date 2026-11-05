# evaluate_v2 verification against frozen e09b (Phase 1d)

Runs: 12 crops x 3 routes x size filter (none, 500) = 72 evaluations.

| check | passing |
|---|---|
| A_markers_dominant_eq_e09b_seed_ignore | 72/72 |
| B_permuted_labels_eq_A | 72/72 |
| C_no_markers_fraction050_eq_e09b_overlap_ignore | 72/72 |
| D_markers_fraction050_eq_e09b_seed_ignore | 72/72 |

A/C/D compare counts, P, R, F1 at IoU 0.25/0.50/0.75 and mean matched IoU at 0.50/0.75 (exact; floats
within 1e-12). mean matched IoU at 0.25 is excluded: e09b keeps an arbitrary maximum-cardinality matching
(label-order dependent); evaluate_v2 takes the max-total-IoU one among them. v2 - e09b at 0.25:
min -0.000000, max +0.011197, nonzero in 62/72 runs.
B compares everything, including mean matched IoU at 0.25.

## Crop-macro F1 (sample SD), by mode

| route | size_filter | mode | f1_iou050_mean | f1_iou050_sd_ddof1 | f1_iou050_median | f1_iou050_min | f1_iou075_mean | f1_iou025_mean | ap_iou050_mean |
|---|---|---|---|---|---|---|---|---|---|
| BA-TOS_ToS | none | markers_dominant | 0.8238 | 0.3670 | 0.9850 | 0.0121 | 0.8175 | 0.8377 | 0.8055 |
| BA-TOS_ToS | none | no_markers_fraction050 | 0.8162 | 0.3616 | 0.9742 | 0.0121 | 0.8091 | 0.8295 | 0.7901 |
| BA-TOS_ToS | none | markers_fraction050 | 0.8238 | 0.3670 | 0.9850 | 0.0121 | 0.8175 | 0.8377 | 0.8055 |
| BA-TOS_ToS | 500 | markers_dominant | 0.8297 | 0.3605 | 0.9869 | 0.0153 | 0.8223 | 0.8471 | 0.8121 |
| BA-TOS_ToS | 500 | no_markers_fraction050 | 0.8225 | 0.3542 | 0.9770 | 0.0153 | 0.8140 | 0.8390 | 0.7969 |
| BA-TOS_ToS | 500 | markers_fraction050 | 0.8297 | 0.3605 | 0.9869 | 0.0153 | 0.8223 | 0.8471 | 0.8121 |
| BA-TOS_sel | none | markers_dominant | 0.9063 | 0.2099 | 0.9850 | 0.2599 | 0.8992 | 0.9375 | 0.8710 |
| BA-TOS_sel | none | no_markers_fraction050 | 0.8947 | 0.2197 | 0.9742 | 0.2130 | 0.8885 | 0.9267 | 0.8536 |
| BA-TOS_sel | none | markers_fraction050 | 0.9063 | 0.2099 | 0.9850 | 0.2599 | 0.8992 | 0.9375 | 0.8710 |
| BA-TOS_sel | 500 | markers_dominant | 0.9088 | 0.2108 | 0.9869 | 0.2599 | 0.9017 | 0.9400 | 0.8758 |
| BA-TOS_sel | 500 | no_markers_fraction050 | 0.8972 | 0.2206 | 0.9770 | 0.2130 | 0.8909 | 0.9293 | 0.8583 |
| BA-TOS_sel | 500 | markers_fraction050 | 0.9088 | 0.2108 | 0.9869 | 0.2599 | 0.9017 | 0.9400 | 0.8758 |
| Otsu+WS(edt_otsu) | none | markers_dominant | 0.8301 | 0.1819 | 0.8816 | 0.2599 | 0.7698 | 0.8959 | 0.7366 |
| Otsu+WS(edt_otsu) | none | no_markers_fraction050 | 0.8120 | 0.1902 | 0.8607 | 0.2130 | 0.7568 | 0.8845 | 0.7113 |
| Otsu+WS(edt_otsu) | none | markers_fraction050 | 0.8301 | 0.1819 | 0.8816 | 0.2599 | 0.7698 | 0.8959 | 0.7366 |
| Otsu+WS(edt_otsu) | 500 | markers_dominant | 0.8301 | 0.1819 | 0.8816 | 0.2599 | 0.7698 | 0.8959 | 0.7366 |
| Otsu+WS(edt_otsu) | 500 | no_markers_fraction050 | 0.8120 | 0.1902 | 0.8607 | 0.2130 | 0.7568 | 0.8845 | 0.7113 |
| Otsu+WS(edt_otsu) | 500 | markers_fraction050 | 0.8301 | 0.1819 | 0.8816 | 0.2599 | 0.7698 | 0.8959 | 0.7366 |

## Pooled object, by mode

| route | size_filter | mode | tp_iou050 | fp_iou050 | fn_iou050 | f1_iou050 | f1_iou075 | ap_iou050 |
|---|---|---|---|---|---|---|---|---|
| BA-TOS_ToS | none | markers_dominant | 2505 | 183 | 303 | 0.9116 | 0.9065 | 0.8375 |
| BA-TOS_ToS | none | no_markers_fraction050 | 2456 | 181 | 352 | 0.9021 | 0.8966 | 0.8217 |
| BA-TOS_ToS | none | markers_fraction050 | 2505 | 183 | 303 | 0.9116 | 0.9065 | 0.8375 |
| BA-TOS_ToS | 500 | markers_dominant | 2505 | 45 | 303 | 0.9351 | 0.9298 | 0.8780 |
| BA-TOS_ToS | 500 | no_markers_fraction050 | 2456 | 43 | 352 | 0.9256 | 0.9199 | 0.8615 |
| BA-TOS_ToS | 500 | markers_fraction050 | 2505 | 45 | 303 | 0.9351 | 0.9298 | 0.8780 |
| BA-TOS_sel | none | markers_dominant | 2630 | 101 | 178 | 0.9496 | 0.9442 | 0.9041 |
| BA-TOS_sel | none | no_markers_fraction050 | 2574 | 94 | 234 | 0.9401 | 0.9350 | 0.8870 |
| BA-TOS_sel | none | markers_fraction050 | 2630 | 101 | 178 | 0.9496 | 0.9442 | 0.9041 |
| BA-TOS_sel | 500 | markers_dominant | 2630 | 86 | 178 | 0.9522 | 0.9468 | 0.9088 |
| BA-TOS_sel | 500 | no_markers_fraction050 | 2574 | 79 | 234 | 0.9427 | 0.9376 | 0.8916 |
| BA-TOS_sel | 500 | markers_fraction050 | 2630 | 86 | 178 | 0.9522 | 0.9468 | 0.9088 |
| Otsu+WS(edt_otsu) | none | markers_dominant | 2547 | 529 | 261 | 0.8657 | 0.7995 | 0.7633 |
| Otsu+WS(edt_otsu) | none | no_markers_fraction050 | 2432 | 493 | 376 | 0.8484 | 0.7874 | 0.7367 |
| Otsu+WS(edt_otsu) | none | markers_fraction050 | 2547 | 529 | 261 | 0.8657 | 0.7995 | 0.7633 |
| Otsu+WS(edt_otsu) | 500 | markers_dominant | 2547 | 529 | 261 | 0.8657 | 0.7995 | 0.7633 |
| Otsu+WS(edt_otsu) | 500 | no_markers_fraction050 | 2432 | 493 | 376 | 0.8484 | 0.7874 | 0.7367 |
| Otsu+WS(edt_otsu) | 500 | markers_fraction050 | 2547 | 529 | 261 | 0.8657 | 0.7995 | 0.7633 |

## Morphology of objects matched at IoU 0.50 (markers_dominant; pooled objects, voxel units)

| route | size_filter | n_matched | volume_rel_error_median | volume_rel_error_median_abs | d_eq_error_median | centroid_error_median | centroid_error_p95 | surface_area_rel_error_median | sphericity_error_median | elongation_error_median | gt_d_eq_median | gt_sphericity_median |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BA-TOS_ToS | none | 2505 | -0.0520 | 0.0596 | -0.5058 | 0.0745 | 8.3918 | -0.0051 | -0.0272 | 0.0054 | 28.6700 | 0.8697 |
| BA-TOS_ToS | 500 | 2505 | -0.0520 | 0.0596 | -0.5058 | 0.0745 | 8.3918 | -0.0051 | -0.0272 | 0.0054 | 28.6700 | 0.8697 |
| BA-TOS_sel | none | 2630 | -0.0524 | 0.0604 | -0.5091 | 0.0771 | 8.5297 | -0.0054 | -0.0271 | 0.0062 | 28.6715 | 0.8695 |
| BA-TOS_sel | 500 | 2630 | -0.0524 | 0.0604 | -0.5091 | 0.0771 | 8.5297 | -0.0054 | -0.0271 | 0.0062 | 28.6715 | 0.8695 |
| Otsu+WS(edt_otsu) | none | 2547 | -0.0536 | 0.0644 | -0.5225 | 0.0875 | 10.4403 | -0.0064 | -0.0285 | 0.0095 | 28.6731 | 0.8695 |
| Otsu+WS(edt_otsu) | 500 | 2547 | -0.0536 | 0.0644 | -0.5225 | 0.0875 | 10.4403 | -0.0064 | -0.0285 | 0.0095 | 28.6731 | 0.8695 |
