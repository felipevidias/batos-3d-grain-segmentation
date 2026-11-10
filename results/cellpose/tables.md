# Phase 2 — Cellpose vs BA-TOS v3, no-marker evaluation (evaluate_v2, fraction 0.50)

v3 = BA-TOS_sel route of `batos-paper-v3-scalefree` (gate never fires on these crops, so sel = ToS).
Crop-macro and pooled tables are separate. Size filter 500 = labels < 500 voxels removed, no relabelling.

## Crop-macro F1, size filter none

| scope | method | n | F1@0.25 | F1@0.50 mean ± SD | median | min (crop) | # < 0.90 | F1@0.75 | mIoU@0.50 | failed |
|---|---|---|---|---|---|---|---|---|---|---|
| dev (12 EFR) | cellpose:cpsam_v2 | 12 | 0.998 | 0.998 ± 0.004 | 1.000 | 0.987 (EFRGP01_01_crop07) | 0 | 0.998 | 0.888 | 0 |
| dev (12 EFR) | cellpose:cpsam | 12 | 1.000 | 1.000 ± 0.001 | 1.000 | 0.998 (EFRGP01_01_crop07) | 0 | 1.000 | 0.891 | 0 |
| dev (12 EFR) | cellpose:cpdino | 12 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.994 (EFRGP01_01_crop07) | 0 | 0.999 | 0.912 | 0 |
| dev (12 EFR) | cellpose:cpdino-vitb | 12 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.995 (EFRGP01_00_crop03) | 0 | 0.999 | 0.896 | 0 |
| dev (12 EFR) | v3 | 12 | 0.985 | 0.979 ± 0.011 | 0.980 | 0.958 (EFRGP01_00_crop06) | 0 | 0.977 | 0.916 | 0 |
| test: all 52 | cellpose:cpsam_v2 | 52 | 0.998 | 0.998 ± 0.004 | 1.000 | 0.980 (EORGP01_00_test01) | 0 | 0.997 | 0.899 | 0 |
| test: all 52 | cellpose:cpsam | 52 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.997 | 0.912 | 0 |
| test: all 52 | cellpose:cpdino | 52 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.928 | 0 |
| test: all 52 | cellpose:cpdino-vitb | 52 | 0.999 | 0.998 ± 0.004 | 1.000 | 0.977 (EMRGP01_00_test02) | 0 | 0.997 | 0.910 | 0 |
| test: all 52 | v3 | 52 | 0.936 | 0.896 ± 0.168 | 0.970 | 0.398 (EORGP01_01_test05) | 10 | 0.853 | 0.902 | 0 |
| test: EFR_internal | cellpose:cpsam_v2 | 12 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.995 (EFRGP01_00_test06) | 0 | 0.992 | 0.888 | 0 |
| test: EFR_internal | cellpose:cpsam | 12 | 0.999 | 0.999 ± 0.001 | 1.000 | 0.996 (EFRGP01_01_test09) | 0 | 0.992 | 0.892 | 0 |
| test: EFR_internal | cellpose:cpdino | 12 | 0.999 | 0.999 ± 0.003 | 1.000 | 0.991 (EFRGP01_01_test10) | 0 | 0.994 | 0.913 | 0 |
| test: EFR_internal | cellpose:cpdino-vitb | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_test01) | 0 | 0.995 | 0.896 | 0 |
| test: EFR_internal | v3 | 12 | 0.982 | 0.977 ± 0.015 | 0.979 | 0.953 (EFRGP01_01_test11) | 0 | 0.966 | 0.909 | 0 |
| test: external | cellpose:cpsam_v2 | 40 | 0.998 | 0.998 ± 0.004 | 1.000 | 0.980 (EORGP01_00_test01) | 0 | 0.998 | 0.903 | 0 |
| test: external | cellpose:cpsam | 40 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.999 | 0.918 | 0 |
| test: external | cellpose:cpdino | 40 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.999 | 0.932 | 0 |
| test: external | cellpose:cpdino-vitb | 40 | 0.998 | 0.998 ± 0.005 | 1.000 | 0.977 (EMRGP01_00_test02) | 0 | 0.998 | 0.915 | 0 |
| test: external | v3 | 40 | 0.922 | 0.872 ± 0.185 | 0.965 | 0.398 (EORGP01_01_test05) | 10 | 0.819 | 0.900 | 0 |
| test: EFR | cellpose:cpsam_v2 | 12 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.995 (EFRGP01_00_test06) | 0 | 0.992 | 0.888 | 0 |
| test: EFR | cellpose:cpsam | 12 | 0.999 | 0.999 ± 0.001 | 1.000 | 0.996 (EFRGP01_01_test09) | 0 | 0.992 | 0.892 | 0 |
| test: EFR | cellpose:cpdino | 12 | 0.999 | 0.999 ± 0.003 | 1.000 | 0.991 (EFRGP01_01_test10) | 0 | 0.994 | 0.913 | 0 |
| test: EFR | cellpose:cpdino-vitb | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_test01) | 0 | 0.995 | 0.896 | 0 |
| test: EFR | v3 | 12 | 0.982 | 0.977 ± 0.015 | 0.979 | 0.953 (EFRGP01_01_test11) | 0 | 0.966 | 0.909 | 0 |
| test: EFS | cellpose:cpsam_v2 | 8 | 0.999 | 0.999 ± 0.003 | 1.000 | 0.993 (EFSGP01_01_test05) | 0 | 0.999 | 0.895 | 0 |
| test: EFS | cellpose:cpsam | 8 | 0.998 | 0.998 ± 0.003 | 1.000 | 0.993 (EFSGP01_01_test05) | 0 | 0.998 | 0.899 | 0 |
| test: EFS | cellpose:cpdino | 8 | 0.999 | 0.999 ± 0.001 | 1.000 | 0.997 (EFSGP01_00_test02) | 0 | 0.999 | 0.917 | 0 |
| test: EFS | cellpose:cpdino-vitb | 8 | 0.999 | 0.999 ± 0.001 | 1.000 | 0.997 (EFSGP01_00_test02) | 0 | 0.999 | 0.903 | 0 |
| test: EFS | v3 | 8 | 0.979 | 0.964 ± 0.048 | 0.986 | 0.858 (EFSGP01_01_test05) | 1 | 0.957 | 0.898 | 0 |
| test: EMR | cellpose:cpsam_v2 | 8 | 0.998 | 0.998 ± 0.003 | 1.000 | 0.991 (EMRGP01_00_test01) | 0 | 0.998 | 0.892 | 0 |
| test: EMR | cellpose:cpsam | 8 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.996 (EMRGP01_01_test08) | 0 | 0.999 | 0.912 | 0 |
| test: EMR | cellpose:cpdino | 8 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.996 (EMRGP01_01_test08) | 0 | 0.999 | 0.928 | 0 |
| test: EMR | cellpose:cpdino-vitb | 8 | 0.996 | 0.996 ± 0.008 | 1.000 | 0.977 (EMRGP01_00_test02) | 0 | 0.996 | 0.909 | 0 |
| test: EMR | v3 | 8 | 0.979 | 0.966 ± 0.032 | 0.978 | 0.910 (EMRGP01_00_test03) | 0 | 0.871 | 0.907 | 0 |
| test: EMS | cellpose:cpsam_v2 | 8 | 0.999 | 0.999 ± 0.003 | 1.000 | 0.994 (EMSGP01_01_test07) | 0 | 0.999 | 0.904 | 0 |
| test: EMS | cellpose:cpsam | 8 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.996 (EMSGP01_01_test08) | 0 | 0.999 | 0.921 | 0 |
| test: EMS | cellpose:cpdino | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EMSGP01_00_test01) | 0 | 1.000 | 0.933 | 0 |
| test: EMS | cellpose:cpdino-vitb | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EMSGP01_00_test01) | 0 | 1.000 | 0.920 | 0 |
| test: EMS | v3 | 8 | 0.978 | 0.970 ± 0.023 | 0.975 | 0.931 (EMSGP01_01_test05) | 0 | 0.960 | 0.925 | 0 |
| test: EOR | cellpose:cpsam_v2 | 8 | 0.997 | 0.996 ± 0.007 | 1.000 | 0.980 (EORGP01_00_test01) | 0 | 0.996 | 0.913 | 0 |
| test: EOR | cellpose:cpsam | 8 | 0.999 | 0.998 ± 0.004 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.932 | 0 |
| test: EOR | cellpose:cpdino | 8 | 0.999 | 0.998 ± 0.004 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.942 | 0 |
| test: EOR | cellpose:cpdino-vitb | 8 | 0.995 | 0.995 ± 0.006 | 0.998 | 0.985 (EORGP01_00_test01) | 0 | 0.995 | 0.922 | 0 |
| test: EOR | v3 | 8 | 0.839 | 0.735 ± 0.243 | 0.748 | 0.398 (EORGP01_01_test05) | 4 | 0.664 | 0.892 | 0 |
| test: EOS | cellpose:cpsam_v2 | 8 | 0.998 | 0.998 ± 0.004 | 1.000 | 0.988 (EOSGP01_00_test03) | 0 | 0.998 | 0.909 | 0 |
| test: EOS | cellpose:cpsam | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EOSGP01_00_test01) | 0 | 1.000 | 0.926 | 0 |
| test: EOS | cellpose:cpdino | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EOSGP01_00_test01) | 0 | 1.000 | 0.939 | 0 |
| test: EOS | cellpose:cpdino-vitb | 8 | 1.000 | 1.000 ± 0.001 | 1.000 | 0.997 (EOSGP01_00_test03) | 0 | 1.000 | 0.921 | 0 |
| test: EOS | v3 | 8 | 0.837 | 0.727 ± 0.226 | 0.727 | 0.421 (EOSGP01_01_test05) | 5 | 0.644 | 0.880 | 0 |

## Pooled, size filter none

| scope | method | GT | pred | TP@0.50 | FP@0.50 | FN@0.50 | P@0.50 | R@0.50 | F1@0.25 | F1@0.50 | F1@0.75 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| dev (12 EFR) | cellpose:cpsam_v2 | 2808 | 2821 | 2808 | 13 | 0 | 0.995 | 1.000 | 0.998 | 0.998 | 0.998 |
| dev (12 EFR) | cellpose:cpsam | 2808 | 2810 | 2808 | 2 | 0 | 0.999 | 1.000 | 1.000 | 1.000 | 1.000 |
| dev (12 EFR) | cellpose:cpdino | 2808 | 2817 | 2808 | 9 | 0 | 0.997 | 1.000 | 0.998 | 0.998 | 0.998 |
| dev (12 EFR) | cellpose:cpdino-vitb | 2808 | 2815 | 2808 | 7 | 0 | 0.998 | 1.000 | 0.999 | 0.999 | 0.999 |
| dev (12 EFR) | v3 | 2808 | 2736 | 2711 | 25 | 97 | 0.991 | 0.965 | 0.984 | 0.978 | 0.975 |
| test: all 52 | cellpose:cpsam_v2 | 8614 | 8645 | 8613 | 32 | 1 | 0.996 | 1.000 | 0.998 | 0.998 | 0.997 |
| test: all 52 | cellpose:cpsam | 8614 | 8630 | 8613 | 17 | 1 | 0.998 | 1.000 | 0.999 | 0.999 | 0.998 |
| test: all 52 | cellpose:cpdino | 8614 | 8625 | 8613 | 12 | 1 | 0.999 | 1.000 | 0.999 | 0.999 | 0.999 |
| test: all 52 | cellpose:cpdino-vitb | 8614 | 8631 | 8613 | 18 | 1 | 0.998 | 1.000 | 0.999 | 0.999 | 0.998 |
| test: all 52 | v3 | 8614 | 7791 | 7437 | 354 | 1177 | 0.955 | 0.863 | 0.943 | 0.907 | 0.875 |
| test: EFR_internal | cellpose:cpsam_v2 | 2476 | 2481 | 2476 | 5 | 0 | 0.998 | 1.000 | 0.999 | 0.999 | 0.996 |
| test: EFR_internal | cellpose:cpsam | 2476 | 2480 | 2476 | 4 | 0 | 0.998 | 1.000 | 0.999 | 0.999 | 0.996 |
| test: EFR_internal | cellpose:cpdino | 2476 | 2480 | 2476 | 4 | 0 | 0.998 | 1.000 | 0.999 | 0.999 | 0.997 |
| test: EFR_internal | cellpose:cpdino-vitb | 2476 | 2476 | 2476 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.998 |
| test: EFR_internal | v3 | 2476 | 2428 | 2394 | 34 | 82 | 0.986 | 0.967 | 0.983 | 0.976 | 0.969 |
| test: external | cellpose:cpsam_v2 | 6138 | 6164 | 6137 | 27 | 1 | 0.996 | 1.000 | 0.998 | 0.998 | 0.998 |
| test: external | cellpose:cpsam | 6138 | 6150 | 6137 | 13 | 1 | 0.998 | 1.000 | 0.999 | 0.999 | 0.999 |
| test: external | cellpose:cpdino | 6138 | 6145 | 6137 | 8 | 1 | 0.999 | 1.000 | 0.999 | 0.999 | 0.999 |
| test: external | cellpose:cpdino-vitb | 6138 | 6155 | 6137 | 18 | 1 | 0.997 | 1.000 | 0.999 | 0.998 | 0.998 |
| test: external | v3 | 6138 | 5363 | 5043 | 320 | 1095 | 0.940 | 0.822 | 0.926 | 0.877 | 0.834 |
| test: EFR | cellpose:cpsam_v2 | 2476 | 2481 | 2476 | 5 | 0 | 0.998 | 1.000 | 0.999 | 0.999 | 0.996 |
| test: EFR | cellpose:cpsam | 2476 | 2480 | 2476 | 4 | 0 | 0.998 | 1.000 | 0.999 | 0.999 | 0.996 |
| test: EFR | cellpose:cpdino | 2476 | 2480 | 2476 | 4 | 0 | 0.998 | 1.000 | 0.999 | 0.999 | 0.997 |
| test: EFR | cellpose:cpdino-vitb | 2476 | 2476 | 2476 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.998 |
| test: EFR | v3 | 2476 | 2428 | 2394 | 34 | 82 | 0.986 | 0.967 | 0.983 | 0.976 | 0.969 |
| test: EFS | cellpose:cpsam_v2 | 1793 | 1798 | 1793 | 5 | 0 | 0.997 | 1.000 | 0.999 | 0.999 | 0.999 |
| test: EFS | cellpose:cpsam | 1793 | 1799 | 1793 | 6 | 0 | 0.997 | 1.000 | 0.998 | 0.998 | 0.998 |
| test: EFS | cellpose:cpdino | 1793 | 1795 | 1793 | 2 | 0 | 0.999 | 1.000 | 0.999 | 0.999 | 0.999 |
| test: EFS | cellpose:cpdino-vitb | 1793 | 1795 | 1793 | 2 | 0 | 0.999 | 1.000 | 0.999 | 0.999 | 0.999 |
| test: EFS | v3 | 1793 | 1720 | 1686 | 34 | 107 | 0.980 | 0.940 | 0.977 | 0.960 | 0.951 |
| test: EMR | cellpose:cpsam_v2 | 1039 | 1044 | 1039 | 5 | 0 | 0.995 | 1.000 | 0.998 | 0.998 | 0.998 |
| test: EMR | cellpose:cpsam | 1039 | 1041 | 1039 | 2 | 0 | 0.998 | 1.000 | 0.999 | 0.999 | 0.999 |
| test: EMR | cellpose:cpdino | 1039 | 1041 | 1039 | 2 | 0 | 0.998 | 1.000 | 0.999 | 0.999 | 0.999 |
| test: EMR | cellpose:cpdino-vitb | 1039 | 1044 | 1039 | 5 | 0 | 0.995 | 1.000 | 0.998 | 0.998 | 0.998 |
| test: EMR | v3 | 1039 | 990 | 973 | 17 | 66 | 0.983 | 0.936 | 0.976 | 0.959 | 0.886 |
| test: EMS | cellpose:cpsam_v2 | 1070 | 1074 | 1070 | 4 | 0 | 0.996 | 1.000 | 0.998 | 0.998 | 0.998 |
| test: EMS | cellpose:cpsam | 1070 | 1071 | 1070 | 1 | 0 | 0.999 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMS | cellpose:cpdino | 1070 | 1070 | 1070 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMS | cellpose:cpdino-vitb | 1070 | 1070 | 1070 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMS | v3 | 1070 | 1018 | 1007 | 11 | 63 | 0.989 | 0.941 | 0.975 | 0.965 | 0.957 |
| test: EOR | cellpose:cpsam_v2 | 1120 | 1127 | 1119 | 8 | 1 | 0.993 | 0.999 | 0.997 | 0.996 | 0.996 |
| test: EOR | cellpose:cpsam | 1120 | 1123 | 1119 | 4 | 1 | 0.996 | 0.999 | 0.999 | 0.998 | 0.998 |
| test: EOR | cellpose:cpdino | 1120 | 1123 | 1119 | 4 | 1 | 0.996 | 0.999 | 0.999 | 0.998 | 0.998 |
| test: EOR | cellpose:cpdino-vitb | 1120 | 1129 | 1119 | 10 | 1 | 0.991 | 0.999 | 0.996 | 0.995 | 0.995 |
| test: EOR | v3 | 1120 | 818 | 692 | 126 | 428 | 0.846 | 0.618 | 0.827 | 0.714 | 0.645 |
| test: EOS | cellpose:cpsam_v2 | 1116 | 1121 | 1116 | 5 | 0 | 0.996 | 1.000 | 0.998 | 0.998 | 0.998 |
| test: EOS | cellpose:cpsam | 1116 | 1116 | 1116 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EOS | cellpose:cpdino | 1116 | 1116 | 1116 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EOS | cellpose:cpdino-vitb | 1116 | 1117 | 1116 | 1 | 0 | 0.999 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EOS | v3 | 1116 | 817 | 685 | 132 | 431 | 0.838 | 0.614 | 0.827 | 0.709 | 0.624 |

## Crop-macro F1, size filter 500

| scope | method | n | F1@0.25 | F1@0.50 mean ± SD | median | min (crop) | # < 0.90 | F1@0.75 | mIoU@0.50 | failed |
|---|---|---|---|---|---|---|---|---|---|---|
| dev (12 EFR) | cellpose:cpsam_v2 | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_crop02) | 0 | 1.000 | 0.888 | 0 |
| dev (12 EFR) | cellpose:cpsam | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_crop02) | 0 | 1.000 | 0.891 | 0 |
| dev (12 EFR) | cellpose:cpdino | 12 | 1.000 | 1.000 ± 0.001 | 1.000 | 0.998 (EFRGP01_00_crop07) | 0 | 1.000 | 0.912 | 0 |
| dev (12 EFR) | cellpose:cpdino-vitb | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_crop02) | 0 | 1.000 | 0.896 | 0 |
| dev (12 EFR) | v3 | 12 | 0.986 | 0.981 ± 0.011 | 0.980 | 0.958 (EFRGP01_00_crop06) | 0 | 0.978 | 0.916 | 0 |
| test: all 52 | cellpose:cpsam_v2 | 52 | 1.000 | 1.000 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.899 | 0 |
| test: all 52 | cellpose:cpsam | 52 | 1.000 | 1.000 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.912 | 0 |
| test: all 52 | cellpose:cpdino | 52 | 1.000 | 1.000 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.928 | 0 |
| test: all 52 | cellpose:cpdino-vitb | 52 | 1.000 | 1.000 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.910 | 0 |
| test: all 52 | v3 | 52 | 0.938 | 0.898 ± 0.168 | 0.972 | 0.398 (EORGP01_01_test05) | 10 | 0.855 | 0.902 | 0 |
| test: EFR_internal | cellpose:cpsam_v2 | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_test01) | 0 | 0.993 | 0.888 | 0 |
| test: EFR_internal | cellpose:cpsam | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_test01) | 0 | 0.992 | 0.892 | 0 |
| test: EFR_internal | cellpose:cpdino | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_test01) | 0 | 0.995 | 0.913 | 0 |
| test: EFR_internal | cellpose:cpdino-vitb | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_test01) | 0 | 0.995 | 0.896 | 0 |
| test: EFR_internal | v3 | 12 | 0.988 | 0.982 ± 0.014 | 0.986 | 0.953 (EFRGP01_01_test11) | 0 | 0.971 | 0.909 | 0 |
| test: external | cellpose:cpsam_v2 | 40 | 1.000 | 1.000 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 1.000 | 0.903 | 0 |
| test: external | cellpose:cpsam | 40 | 1.000 | 1.000 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 1.000 | 0.918 | 0 |
| test: external | cellpose:cpdino | 40 | 1.000 | 1.000 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 1.000 | 0.932 | 0 |
| test: external | cellpose:cpdino-vitb | 40 | 1.000 | 1.000 ± 0.002 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 1.000 | 0.915 | 0 |
| test: external | v3 | 40 | 0.923 | 0.873 ± 0.184 | 0.965 | 0.398 (EORGP01_01_test05) | 10 | 0.820 | 0.900 | 0 |
| test: EFR | cellpose:cpsam_v2 | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_test01) | 0 | 0.993 | 0.888 | 0 |
| test: EFR | cellpose:cpsam | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_test01) | 0 | 0.992 | 0.892 | 0 |
| test: EFR | cellpose:cpdino | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_test01) | 0 | 0.995 | 0.913 | 0 |
| test: EFR | cellpose:cpdino-vitb | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_test01) | 0 | 0.995 | 0.896 | 0 |
| test: EFR | v3 | 12 | 0.988 | 0.982 ± 0.014 | 0.986 | 0.953 (EFRGP01_01_test11) | 0 | 0.971 | 0.909 | 0 |
| test: EFS | cellpose:cpsam_v2 | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFSGP01_00_test01) | 0 | 1.000 | 0.895 | 0 |
| test: EFS | cellpose:cpsam | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFSGP01_00_test01) | 0 | 1.000 | 0.899 | 0 |
| test: EFS | cellpose:cpdino | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFSGP01_00_test01) | 0 | 1.000 | 0.917 | 0 |
| test: EFS | cellpose:cpdino-vitb | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFSGP01_00_test01) | 0 | 1.000 | 0.903 | 0 |
| test: EFS | v3 | 8 | 0.980 | 0.965 ± 0.048 | 0.987 | 0.860 (EFSGP01_01_test05) | 1 | 0.958 | 0.898 | 0 |
| test: EMR | cellpose:cpsam_v2 | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EMRGP01_00_test01) | 0 | 1.000 | 0.892 | 0 |
| test: EMR | cellpose:cpsam | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EMRGP01_00_test01) | 0 | 1.000 | 0.912 | 0 |
| test: EMR | cellpose:cpdino | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EMRGP01_00_test01) | 0 | 1.000 | 0.928 | 0 |
| test: EMR | cellpose:cpdino-vitb | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EMRGP01_00_test01) | 0 | 1.000 | 0.909 | 0 |
| test: EMR | v3 | 8 | 0.979 | 0.966 ± 0.032 | 0.978 | 0.910 (EMRGP01_00_test03) | 0 | 0.871 | 0.907 | 0 |
| test: EMS | cellpose:cpsam_v2 | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EMSGP01_00_test01) | 0 | 1.000 | 0.904 | 0 |
| test: EMS | cellpose:cpsam | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EMSGP01_00_test01) | 0 | 1.000 | 0.921 | 0 |
| test: EMS | cellpose:cpdino | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EMSGP01_00_test01) | 0 | 1.000 | 0.933 | 0 |
| test: EMS | cellpose:cpdino-vitb | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EMSGP01_00_test01) | 0 | 1.000 | 0.920 | 0 |
| test: EMS | v3 | 8 | 0.978 | 0.970 ± 0.023 | 0.975 | 0.931 (EMSGP01_01_test05) | 0 | 0.960 | 0.925 | 0 |
| test: EOR | cellpose:cpsam_v2 | 8 | 0.999 | 0.998 ± 0.004 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.913 | 0 |
| test: EOR | cellpose:cpsam | 8 | 0.999 | 0.998 ± 0.004 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.932 | 0 |
| test: EOR | cellpose:cpdino | 8 | 0.999 | 0.998 ± 0.004 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.942 | 0 |
| test: EOR | cellpose:cpdino-vitb | 8 | 0.999 | 0.998 ± 0.004 | 1.000 | 0.988 (EORGP01_00_test01) | 0 | 0.998 | 0.922 | 0 |
| test: EOR | v3 | 8 | 0.839 | 0.736 ± 0.243 | 0.748 | 0.398 (EORGP01_01_test05) | 4 | 0.665 | 0.892 | 0 |
| test: EOS | cellpose:cpsam_v2 | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EOSGP01_00_test01) | 0 | 1.000 | 0.909 | 0 |
| test: EOS | cellpose:cpsam | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EOSGP01_00_test01) | 0 | 1.000 | 0.926 | 0 |
| test: EOS | cellpose:cpdino | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EOSGP01_00_test01) | 0 | 1.000 | 0.939 | 0 |
| test: EOS | cellpose:cpdino-vitb | 8 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EOSGP01_00_test01) | 0 | 1.000 | 0.921 | 0 |
| test: EOS | v3 | 8 | 0.840 | 0.729 ± 0.225 | 0.729 | 0.423 (EOSGP01_01_test05) | 5 | 0.646 | 0.880 | 0 |

## Pooled, size filter 500

| scope | method | GT | pred | TP@0.50 | FP@0.50 | FN@0.50 | P@0.50 | R@0.50 | F1@0.25 | F1@0.50 | F1@0.75 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| dev (12 EFR) | cellpose:cpsam_v2 | 2808 | 2808 | 2808 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| dev (12 EFR) | cellpose:cpsam | 2808 | 2808 | 2808 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| dev (12 EFR) | cellpose:cpdino | 2808 | 2809 | 2808 | 1 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| dev (12 EFR) | cellpose:cpdino-vitb | 2808 | 2808 | 2808 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| dev (12 EFR) | v3 | 2808 | 2730 | 2711 | 19 | 97 | 0.993 | 0.965 | 0.985 | 0.979 | 0.976 |
| test: all 52 | cellpose:cpsam_v2 | 8614 | 8617 | 8613 | 4 | 1 | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 |
| test: all 52 | cellpose:cpsam | 8614 | 8617 | 8613 | 4 | 1 | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 |
| test: all 52 | cellpose:cpdino | 8614 | 8617 | 8613 | 4 | 1 | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 |
| test: all 52 | cellpose:cpdino-vitb | 8614 | 8617 | 8613 | 4 | 1 | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 |
| test: all 52 | v3 | 8614 | 7762 | 7437 | 325 | 1177 | 0.958 | 0.863 | 0.945 | 0.908 | 0.876 |
| test: EFR_internal | cellpose:cpsam_v2 | 2476 | 2476 | 2476 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.997 |
| test: EFR_internal | cellpose:cpsam | 2476 | 2476 | 2476 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.996 |
| test: EFR_internal | cellpose:cpdino | 2476 | 2476 | 2476 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.998 |
| test: EFR_internal | cellpose:cpdino-vitb | 2476 | 2476 | 2476 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.998 |
| test: EFR_internal | v3 | 2476 | 2411 | 2394 | 17 | 82 | 0.993 | 0.967 | 0.987 | 0.980 | 0.973 |
| test: external | cellpose:cpsam_v2 | 6138 | 6141 | 6137 | 4 | 1 | 0.999 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: external | cellpose:cpsam | 6138 | 6141 | 6137 | 4 | 1 | 0.999 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: external | cellpose:cpdino | 6138 | 6141 | 6137 | 4 | 1 | 0.999 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: external | cellpose:cpdino-vitb | 6138 | 6141 | 6137 | 4 | 1 | 0.999 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: external | v3 | 6138 | 5351 | 5043 | 308 | 1095 | 0.942 | 0.822 | 0.927 | 0.878 | 0.835 |
| test: EFR | cellpose:cpsam_v2 | 2476 | 2476 | 2476 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.997 |
| test: EFR | cellpose:cpsam | 2476 | 2476 | 2476 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.996 |
| test: EFR | cellpose:cpdino | 2476 | 2476 | 2476 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.998 |
| test: EFR | cellpose:cpdino-vitb | 2476 | 2476 | 2476 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.998 |
| test: EFR | v3 | 2476 | 2411 | 2394 | 17 | 82 | 0.993 | 0.967 | 0.987 | 0.980 | 0.973 |
| test: EFS | cellpose:cpsam_v2 | 1793 | 1793 | 1793 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EFS | cellpose:cpsam | 1793 | 1793 | 1793 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EFS | cellpose:cpdino | 1793 | 1793 | 1793 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EFS | cellpose:cpdino-vitb | 1793 | 1793 | 1793 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EFS | v3 | 1793 | 1717 | 1686 | 31 | 107 | 0.982 | 0.940 | 0.978 | 0.961 | 0.952 |
| test: EMR | cellpose:cpsam_v2 | 1039 | 1039 | 1039 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMR | cellpose:cpsam | 1039 | 1039 | 1039 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMR | cellpose:cpdino | 1039 | 1039 | 1039 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMR | cellpose:cpdino-vitb | 1039 | 1039 | 1039 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMR | v3 | 1039 | 990 | 973 | 17 | 66 | 0.983 | 0.936 | 0.976 | 0.959 | 0.886 |
| test: EMS | cellpose:cpsam_v2 | 1070 | 1070 | 1070 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMS | cellpose:cpsam | 1070 | 1070 | 1070 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMS | cellpose:cpdino | 1070 | 1070 | 1070 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMS | cellpose:cpdino-vitb | 1070 | 1070 | 1070 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EMS | v3 | 1070 | 1018 | 1007 | 11 | 63 | 0.989 | 0.941 | 0.975 | 0.965 | 0.957 |
| test: EOR | cellpose:cpsam_v2 | 1120 | 1123 | 1119 | 4 | 1 | 0.996 | 0.999 | 0.999 | 0.998 | 0.998 |
| test: EOR | cellpose:cpsam | 1120 | 1123 | 1119 | 4 | 1 | 0.996 | 0.999 | 0.999 | 0.998 | 0.998 |
| test: EOR | cellpose:cpdino | 1120 | 1123 | 1119 | 4 | 1 | 0.996 | 0.999 | 0.999 | 0.998 | 0.998 |
| test: EOR | cellpose:cpdino-vitb | 1120 | 1123 | 1119 | 4 | 1 | 0.996 | 0.999 | 0.999 | 0.998 | 0.998 |
| test: EOR | v3 | 1120 | 816 | 692 | 124 | 428 | 0.848 | 0.618 | 0.827 | 0.715 | 0.646 |
| test: EOS | cellpose:cpsam_v2 | 1116 | 1116 | 1116 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EOS | cellpose:cpsam | 1116 | 1116 | 1116 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EOS | cellpose:cpdino | 1116 | 1116 | 1116 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EOS | cellpose:cpdino-vitb | 1116 | 1116 | 1116 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| test: EOS | v3 | 1116 | 810 | 685 | 125 | 431 | 0.846 | 0.614 | 0.830 | 0.711 | 0.626 |

## Runtime and memory per crop (200³)

| split | method | device | n | seconds median [min–max] | peak VRAM reserved MB (max) | peak host MB (max) |
|---|---|---|---|---|---|---|
| dev | cellpose:cpsam_v2 | cuda (NVIDIA GeForce RTX 5060) | 12 | 30 [29–46] | 1814 | 3484 |
| dev | cellpose:cpsam | cuda (NVIDIA GeForce RTX 5060) | 12 | 31 [30–47] | 1814 | 3410 |
| dev | cellpose:cpdino | cuda (NVIDIA GeForce RTX 5060) | 12 | 58 [58–77] | 1652 | 3293 |
| dev | cellpose:cpdino-vitb | cuda (NVIDIA GeForce RTX 5060) | 12 | 23 [23–28] | 1188 | 2482 |
| dev | v3 | cpu (AMD Ryzen 7 5700X) | 12 | 397 [361–484] | — | 3355 |
| test | cellpose:cpsam_v2 | cuda (NVIDIA GeForce RTX 5060) | 52 | 30 [29–32] | 1816 | 3540 |
| test | cellpose:cpsam | cuda (NVIDIA GeForce RTX 5060) | 52 | 31 [30–60] | 1816 | 3473 |
| test | cellpose:cpdino | cuda (NVIDIA GeForce RTX 5060) | 52 | 58 [57–58] | 1654 | 3325 |
| test | cellpose:cpdino-vitb | cuda (NVIDIA GeForce RTX 5060) | 52 | 23 [22–23] | 1192 | 2482 |
| test | v3 | cpu (AMD Ryzen 7 5700X) | 52 | 397 [318–567] | — | 3648 |

Cellpose seconds = model.eval only (model load excluded), GPU. v3 seconds = whole pipeline, CPU, not contention-controlled (from v3_run.json).
