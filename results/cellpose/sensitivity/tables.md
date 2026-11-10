# Phase 2 — Cellpose vs BA-TOS v3, no-marker evaluation (evaluate_v2, fraction 0.50)

v3 = BA-TOS_sel route of `batos-paper-v3-scalefree` (gate never fires on these crops, so sel = ToS).
Crop-macro and pooled tables are separate. Size filter 500 = labels < 500 voxels removed, no relabelling.

## Crop-macro F1, size filter none

| scope | method | n | F1@0.25 | F1@0.50 mean ± SD | median | min (crop) | # < 0.90 | F1@0.75 | mIoU@0.50 | failed |
|---|---|---|---|---|---|---|---|---|---|---|
| dev (12 EFR) | cellpose:cpsam_v2 | 12 | 0.998 | 0.998 ± 0.004 | 1.000 | 0.987 (EFRGP01_01_crop07) | 0 | 0.998 | 0.888 | 0 |
| dev (12 EFR) | cellpose:cpsam_v2__flow3D_smooth-1 | 12 | 0.998 | 0.998 ± 0.003 | 1.000 | 0.989 (EFRGP01_01_crop07) | 0 | 0.998 | 0.888 | 0 |
| dev (12 EFR) | cellpose:cpsam_v2__flow3D_smooth-2 | 12 | 0.999 | 0.999 ± 0.002 | 1.000 | 0.994 (EFRGP01_01_crop07) | 0 | 0.999 | 0.888 | 0 |
| dev (12 EFR) | cellpose:cpsam_v2__diameter-2r_granulo | 12 | 0.997 | 0.997 ± 0.004 | 0.998 | 0.987 (EFRGP01_01_crop06) | 0 | 0.997 | 0.861 | 0 |
| dev (12 EFR) | v3 | 12 | 0.985 | 0.979 ± 0.011 | 0.980 | 0.958 (EFRGP01_00_crop06) | 0 | 0.977 | 0.916 | 0 |

## Pooled, size filter none

| scope | method | GT | pred | TP@0.50 | FP@0.50 | FN@0.50 | P@0.50 | R@0.50 | F1@0.25 | F1@0.50 | F1@0.75 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| dev (12 EFR) | cellpose:cpsam_v2 | 2808 | 2821 | 2808 | 13 | 0 | 0.995 | 1.000 | 0.998 | 0.998 | 0.998 |
| dev (12 EFR) | cellpose:cpsam_v2__flow3D_smooth-1 | 2808 | 2818 | 2808 | 10 | 0 | 0.996 | 1.000 | 0.998 | 0.998 | 0.998 |
| dev (12 EFR) | cellpose:cpsam_v2__flow3D_smooth-2 | 2808 | 2815 | 2808 | 7 | 0 | 0.998 | 1.000 | 0.999 | 0.999 | 0.999 |
| dev (12 EFR) | cellpose:cpsam_v2__diameter-2r_granulo | 2808 | 2827 | 2808 | 19 | 0 | 0.993 | 1.000 | 0.997 | 0.997 | 0.997 |
| dev (12 EFR) | v3 | 2808 | 2736 | 2711 | 25 | 97 | 0.991 | 0.965 | 0.984 | 0.978 | 0.975 |

## Crop-macro F1, size filter 500

| scope | method | n | F1@0.25 | F1@0.50 mean ± SD | median | min (crop) | # < 0.90 | F1@0.75 | mIoU@0.50 | failed |
|---|---|---|---|---|---|---|---|---|---|---|
| dev (12 EFR) | cellpose:cpsam_v2 | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_crop02) | 0 | 1.000 | 0.888 | 0 |
| dev (12 EFR) | cellpose:cpsam_v2__flow3D_smooth-1 | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_crop02) | 0 | 1.000 | 0.888 | 0 |
| dev (12 EFR) | cellpose:cpsam_v2__flow3D_smooth-2 | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_crop02) | 0 | 1.000 | 0.888 | 0 |
| dev (12 EFR) | cellpose:cpsam_v2__diameter-2r_granulo | 12 | 1.000 | 1.000 ± 0.000 | 1.000 | 1.000 (EFRGP01_00_crop02) | 0 | 1.000 | 0.861 | 0 |
| dev (12 EFR) | v3 | 12 | 0.986 | 0.981 ± 0.011 | 0.980 | 0.958 (EFRGP01_00_crop06) | 0 | 0.978 | 0.916 | 0 |

## Pooled, size filter 500

| scope | method | GT | pred | TP@0.50 | FP@0.50 | FN@0.50 | P@0.50 | R@0.50 | F1@0.25 | F1@0.50 | F1@0.75 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| dev (12 EFR) | cellpose:cpsam_v2 | 2808 | 2808 | 2808 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| dev (12 EFR) | cellpose:cpsam_v2__flow3D_smooth-1 | 2808 | 2808 | 2808 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| dev (12 EFR) | cellpose:cpsam_v2__flow3D_smooth-2 | 2808 | 2808 | 2808 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| dev (12 EFR) | cellpose:cpsam_v2__diameter-2r_granulo | 2808 | 2808 | 2808 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| dev (12 EFR) | v3 | 2808 | 2730 | 2711 | 19 | 97 | 0.993 | 0.965 | 0.985 | 0.979 | 0.976 |

## Runtime and memory per crop (200³)

| split | method | device | n | seconds median [min–max] | peak VRAM reserved MB (max) | peak host MB (max) |
|---|---|---|---|---|---|---|
| dev | cellpose:cpsam_v2 | cuda (NVIDIA GeForce RTX 5060) | 12 | 30 [29–46] | 1814 | 3484 |
| dev | cellpose:cpsam_v2__flow3D_smooth-1 | cuda (NVIDIA GeForce RTX 5060) | 12 | 30 [30–32] | 1814 | 3658 |
| dev | cellpose:cpsam_v2__flow3D_smooth-2 | cuda (NVIDIA GeForce RTX 5060) | 12 | 31 [30–32] | 1814 | 3708 |
| dev | cellpose:cpsam_v2__diameter-2r_granulo | cuda (NVIDIA GeForce RTX 5060) | 12 | 171 [170–173] | 1814 | 4076 |
| dev | v3 | cpu (AMD Ryzen 7 5700X) | 12 | 397 [361–484] | — | 3355 |

Cellpose seconds = model.eval only (model load excluded), GPU. v3 seconds = whole pipeline, CPU, not contention-controlled (from v3_run.json).
