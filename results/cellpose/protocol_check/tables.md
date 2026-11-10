# Protocol check: no-marker rule vs STRICT (nothing ignored), size filter off

| scope | method | mode | n | F1@0.50 mean ± SD | min | pooled TP | FP | FN | pooled F1@0.50 | F1@0.75 pooled | pred total | ignored | pred/GT |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dev (12 EFR) | cellpose:cpsam_v2 | no_markers_fraction050 | 12 | 0.998 ± 0.004 | 0.987 | 2808 | 13 | 0 | 0.998 | 0.998 | 6331 | 3510 | 1.005 |
| dev (12 EFR) | cellpose:cpsam_v2 | strict_noignore | 12 | 0.616 ± 0.027 | 0.585 | 2808 | 3519 | 0 | 0.615 | 0.615 | 6331 | 0 | 2.253 |
| dev (12 EFR) | v3 | no_markers_fraction050 | 12 | 0.979 ± 0.011 | 0.958 | 2711 | 25 | 97 | 0.978 | 0.975 | 4078 | 1342 | 0.974 |
| dev (12 EFR) | v3 | strict_noignore | 12 | 0.831 ± 0.019 | 0.802 | 2773 | 1107 | 35 | 0.829 | 0.826 | 4078 | 0 | 1.382 |
| test: all 52 | cellpose:cpsam_v2 | no_markers_fraction050 | 52 | 0.998 ± 0.004 | 0.980 | 8613 | 32 | 1 | 0.998 | 0.997 | 20838 | 12193 | 1.004 |
| test: all 52 | cellpose:cpsam_v2 | strict_noignore | 52 | 0.585 ± 0.031 | 0.538 | 8613 | 12193 | 1 | 0.586 | 0.585 | 20838 | 0 | 2.415 |
| test: all 52 | v3 | no_markers_fraction050 | 52 | 0.896 ± 0.168 | 0.398 | 7437 | 354 | 1177 | 0.907 | 0.875 | 11755 | 3964 | 0.904 |
| test: all 52 | v3 | strict_noignore | 52 | 0.770 ± 0.119 | 0.418 | 7697 | 3516 | 917 | 0.776 | 0.748 | 11755 | 0 | 1.302 |
| test: EFR_internal | cellpose:cpsam_v2 | no_markers_fraction050 | 12 | 0.999 ± 0.002 | 0.995 | 2476 | 5 | 0 | 0.999 | 0.996 | 5756 | 3275 | 1.002 |
| test: EFR_internal | cellpose:cpsam_v2 | strict_noignore | 12 | 0.605 ± 0.023 | 0.576 | 2476 | 3274 | 0 | 0.602 | 0.600 | 5756 | 0 | 2.322 |
| test: EFR_internal | v3 | no_markers_fraction050 | 12 | 0.977 ± 0.015 | 0.953 | 2394 | 34 | 82 | 0.976 | 0.969 | 3614 | 1186 | 0.981 |
| test: EFR_internal | v3 | strict_noignore | 12 | 0.823 ± 0.027 | 0.793 | 2440 | 1025 | 36 | 0.821 | 0.815 | 3614 | 0 | 1.399 |
| test: external | cellpose:cpsam_v2 | no_markers_fraction050 | 40 | 0.998 ± 0.004 | 0.980 | 6137 | 27 | 1 | 0.998 | 0.998 | 15082 | 8918 | 1.004 |
| test: external | cellpose:cpsam_v2 | strict_noignore | 40 | 0.579 ± 0.032 | 0.538 | 6137 | 8919 | 1 | 0.579 | 0.579 | 15082 | 0 | 2.453 |
| test: external | v3 | no_markers_fraction050 | 40 | 0.872 ± 0.185 | 0.398 | 5043 | 320 | 1095 | 0.877 | 0.834 | 8141 | 2778 | 0.874 |
| test: external | v3 | strict_noignore | 40 | 0.754 ± 0.132 | 0.418 | 5257 | 2491 | 881 | 0.757 | 0.720 | 8141 | 0 | 1.262 |
| test: EFR | cellpose:cpsam_v2 | no_markers_fraction050 | 12 | 0.999 ± 0.002 | 0.995 | 2476 | 5 | 0 | 0.999 | 0.996 | 5756 | 3275 | 1.002 |
| test: EFR | cellpose:cpsam_v2 | strict_noignore | 12 | 0.605 ± 0.023 | 0.576 | 2476 | 3274 | 0 | 0.602 | 0.600 | 5756 | 0 | 2.322 |
| test: EFR | v3 | no_markers_fraction050 | 12 | 0.977 ± 0.015 | 0.953 | 2394 | 34 | 82 | 0.976 | 0.969 | 3614 | 1186 | 0.981 |
| test: EFR | v3 | strict_noignore | 12 | 0.823 ± 0.027 | 0.793 | 2440 | 1025 | 36 | 0.821 | 0.815 | 3614 | 0 | 1.399 |
| test: EFS | cellpose:cpsam_v2 | no_markers_fraction050 | 8 | 0.999 ± 0.003 | 0.993 | 1793 | 5 | 0 | 0.999 | 0.999 | 4093 | 2295 | 1.003 |
| test: EFS | cellpose:cpsam_v2 | strict_noignore | 8 | 0.614 ± 0.021 | 0.586 | 1793 | 2294 | 0 | 0.610 | 0.610 | 4093 | 0 | 2.279 |
| test: EFS | v3 | no_markers_fraction050 | 8 | 0.964 ± 0.048 | 0.858 | 1686 | 34 | 107 | 0.960 | 0.951 | 2535 | 815 | 0.959 |
| test: EFS | v3 | strict_noignore | 8 | 0.816 ± 0.036 | 0.748 | 1721 | 727 | 72 | 0.812 | 0.804 | 2535 | 0 | 1.365 |
| test: EMR | cellpose:cpsam_v2 | no_markers_fraction050 | 8 | 0.998 ± 0.003 | 0.991 | 1039 | 5 | 0 | 0.998 | 0.998 | 2647 | 1603 | 1.005 |
| test: EMR | cellpose:cpsam_v2 | strict_noignore | 8 | 0.565 ± 0.025 | 0.545 | 1039 | 1607 | 0 | 0.564 | 0.564 | 2647 | 0 | 2.547 |
| test: EMR | v3 | no_markers_fraction050 | 8 | 0.966 ± 0.032 | 0.910 | 973 | 17 | 66 | 0.959 | 0.886 | 1506 | 516 | 0.953 |
| test: EMR | v3 | strict_noignore | 8 | 0.812 ± 0.038 | 0.768 | 1002 | 447 | 37 | 0.805 | 0.746 | 1506 | 0 | 1.395 |
| test: EMS | cellpose:cpsam_v2 | no_markers_fraction050 | 8 | 0.999 ± 0.003 | 0.994 | 1070 | 4 | 0 | 0.998 | 0.998 | 2761 | 1687 | 1.004 |
| test: EMS | cellpose:cpsam_v2 | strict_noignore | 8 | 0.560 ± 0.024 | 0.538 | 1070 | 1690 | 0 | 0.559 | 0.559 | 2761 | 0 | 2.579 |
| test: EMS | v3 | no_markers_fraction050 | 8 | 0.970 ± 0.023 | 0.931 | 1007 | 11 | 63 | 0.965 | 0.957 | 1619 | 601 | 0.951 |
| test: EMS | v3 | strict_noignore | 8 | 0.800 ± 0.031 | 0.755 | 1046 | 499 | 24 | 0.800 | 0.793 | 1619 | 0 | 1.444 |
| test: EOR | cellpose:cpsam_v2 | no_markers_fraction050 | 8 | 0.996 ± 0.007 | 0.980 | 1119 | 8 | 1 | 0.996 | 0.996 | 2804 | 1677 | 1.006 |
| test: EOR | cellpose:cpsam_v2 | strict_noignore | 8 | 0.576 ± 0.027 | 0.542 | 1119 | 1676 | 1 | 0.572 | 0.572 | 2804 | 0 | 2.496 |
| test: EOR | v3 | no_markers_fraction050 | 8 | 0.735 ± 0.243 | 0.398 | 692 | 126 | 428 | 0.714 | 0.645 | 1219 | 401 | 0.730 |
| test: EOR | v3 | strict_noignore | 8 | 0.674 ± 0.184 | 0.418 | 747 | 399 | 373 | 0.659 | 0.596 | 1219 | 0 | 1.023 |
| test: EOS | cellpose:cpsam_v2 | no_markers_fraction050 | 8 | 0.998 ± 0.004 | 0.988 | 1116 | 5 | 0 | 0.998 | 0.998 | 2777 | 1656 | 1.004 |
| test: EOS | cellpose:cpsam_v2 | strict_noignore | 8 | 0.581 ± 0.034 | 0.538 | 1116 | 1652 | 0 | 0.575 | 0.575 | 2777 | 0 | 2.480 |
| test: EOS | v3 | no_markers_fraction050 | 8 | 0.727 ± 0.226 | 0.421 | 685 | 132 | 431 | 0.709 | 0.624 | 1262 | 445 | 0.732 |
| test: EOS | v3 | strict_noignore | 8 | 0.667 ± 0.180 | 0.435 | 741 | 419 | 375 | 0.651 | 0.574 | 1262 | 0 | 1.039 |
