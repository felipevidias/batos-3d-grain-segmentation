# Phase 1f — test protocol: frozen vs fixflat on 52 held-out crops (6 specimens)

Date 2026-11-06/07. Branch `journal-extension`. Method under test: **batos-paper-v2-fixflat** (tag `batos-paper-v2-fixflat`,
`configs/batos_paper_v2_fixflat.yaml`), compared with **batos-paper-v2-frozen** (with the E23c gate). **No parameter was
changed for or after this test.** GT was used only for evaluation and for the diagnostic columns (GT inscribed radius).

## 1. Protocol (as pre-registered)

* Crops: `data/metadata/crops_test.json` (committed in 4a72daa **before** cutting; GT-free rule, seed 20261106). 52 crops of
  200³, zyx: EFR 12 (6 per scan, no overlap with the 12 dev crops), EFS/EMR/EMS/EOR/EOS 8 each (4 per scan); every
  specimen and scan has ≥ 1 boundary crop (specimen fraction 0.40–0.80 from the greyscale cylinder fit). Manifest with
  SHA-256: `data/metadata/crops_manifest_test.csv` (8,614 valid GT grains).
* Runs: `scripts/run_test_protocol.py`, sequential, peak ≤ 3.7 GB (budget 8 GB). Frozen = `python -I -m batos_v2.frozen`
  (verbatim E-series scripts); fixflat = frozen scale + openrec, then `batos_v2.flatten_fixed --variant A`, then frozen code.
* Evaluation: `batos_v2.evaluate_v2`, IoU 0.25/0.50/0.75, size filter off and 500; modes `markers_dominant` (= e09b
  seed-ignore, paper protocol, **primary**) and `no_markers_fraction050` (own-voxel support, overlap-ignore 0.50).
* Aggregation: `scripts/aggregate_test_protocol.py`; forensics: `scripts/forensics_test_protocol.py`.

### Deviations and implementation notes (all disclosed, none tuned)

1. **e13 needs exactly 12 manifest rows.** Each specimen's crops were run through e13 verbatim, padded to 12 rows with
   duplicate rows (`__pad`); the script checks that every pad reproduces its original row exactly (e13 is per-crop
   independent). Passed for all 6 specimens.
2. **openrec implementation:** `--openrec-impl auto` would switch to the bit-identical `openrec_lowmem` above radius 15. It
   never triggered: **all 51 runs used the frozen `openrec_v2`** (radii 10–15).
3. **Pipeline crashes are method failures, scored as empty segmentations** (F1 = 0, every valid GT grain a FN; no retry).
   Five crashes, all `e11_build_markers_ablation.py`: "Nenhum candidato MIN caiu dentro de internal_dark" (no MIN
   candidate inside internal_dark): fixflat on EMS 00_test01, EOR 00_test01, EOS 00_test01, EOS 01_test05; frozen on
   EOS 00_test03. On EOS 00_test03 fixflat could not run at all (the harness reuses the frozen scale/openrec outputs and
   requires a complete frozen run); it is also scored as a failure (conservative against fixflat). This rule was added to
   the driver when the first crash stopped the run (crop 29); the crop sets and parameters were not touched.
4. EFR GT came from the existing `C:\BA_TOS_DATA\gt_full`; the other 5 specimens' `_00/_01` label volumes were extracted
   with py7zr into `gt_full_tmp`, cut, and deleted.

## 2. Results — crop-macro F1 (markers/dominant = seed-ignore)

`EFR_internal` = 12 new EFR crops; `external` = 40 crops of the 5 other specimens. BA-TOS_ToS = fixflat route without
gate; BA-TOS_sel = with E23c gate. For fixflat ToS ≡ sel except where the gate fires (it never does on EFR/EFS).

### Size filter off

| scope | method | route | n | F1@0.50 mean ± SD | median | min (crop) | F1@0.75 | F1@0.25 |
|---|---|---|---|---|---|---|---|---|
| all | frozen | ToS | 52 | 0.343 ± 0.426 | 0.060 | 0.000 | 0.319 | 0.384 |
| all | frozen | sel | 52 | 0.538 ± 0.437 | 0.735 | 0.000 | 0.492 | 0.576 |
| all | fixflat | ToS | 52 | 0.570 ± 0.470 | 0.919 | 0.000 | 0.542 | 0.582 |
| all | fixflat | sel | 52 | 0.575 ± 0.465 | 0.919 | 0.000 | 0.542 | 0.587 |
| EFR_internal | frozen | ToS | 12 | 0.668 ± 0.465 | 0.973 | 0.015 (00_test06) | 0.661 | 0.686 |
| EFR_internal | frozen | sel | 12 | 0.879 ± 0.223 | 0.973 | 0.212 (00_test02) | 0.846 | 0.915 |
| EFR_internal | fixflat | ToS = sel | 12 | **0.956 ± 0.048** | 0.973 | 0.835 (00_test06) | 0.919 | 0.966 |
| external | frozen | ToS | 40 | 0.245 ± 0.365 | 0.045 | 0.000 | 0.216 | 0.294 |
| external | frozen | sel | 40 | 0.435 ± 0.435 | 0.153 | 0.000 | 0.385 | 0.474 |
| external | fixflat | ToS | 40 | 0.454 ± 0.478 | 0.078 | 0.000 | 0.428 | 0.467 |
| external | fixflat | sel | 40 | 0.460 ± 0.473 | 0.155 | 0.000 | 0.429 | 0.474 |

Per specimen (sel route, F1@0.50 mean ± SD / min; ToS in `tables.md`):

| specimen | frozen sel | fixflat sel | fixflat F1@0.75 |
|---|---|---|---|
| EFR (12) | 0.879 ± 0.223 / 0.212 | **0.956 ± 0.048 / 0.835** | 0.919 |
| EFS (8) | 0.882 ± 0.155 / 0.609 | **0.968 ± 0.043 / 0.874** | 0.961 |
| EMR (8) | 0.303 ± 0.439 / 0.000 | 0.250 ± 0.416 / 0.000 | 0.232 |
| EMS (8) | 0.432 ± 0.465 / 0.000 | 0.479 ± 0.474 / 0.000 | 0.444 |
| EOR (8) | 0.293 ± 0.362 / 0.000 | 0.238 ± 0.434 / 0.000 | 0.205 |
| EOS (8) | 0.265 ± 0.429 / 0.000 | 0.365 ± 0.499 / 0.000 | 0.304 |

### Size filter 500

| scope | method | route | F1@0.50 mean ± SD | median | min | F1@0.75 |
|---|---|---|---|---|---|---|
| all | frozen | sel | 0.555 ± 0.447 | 0.816 | 0.000 | 0.507 |
| all | fixflat | sel | 0.587 ± 0.475 | 0.950 | 0.000 | 0.553 |
| EFR_internal | frozen | sel | 0.883 ± 0.223 | 0.979 | 0.212 | 0.849 |
| EFR_internal | fixflat | ToS = sel | **0.985 ± 0.021** | 0.991 | 0.924 (00_test02) | 0.945 |
| external | frozen | sel | 0.456 ± 0.452 | 0.153 | 0.000 | 0.405 |
| external | fixflat | sel | 0.467 ± 0.481 | 0.155 | 0.000 | 0.436 |

### No-marker mode (own-voxel support, overlap-ignore 0.50), crop-macro F1@0.50, filter off / 500

| scope | frozen sel | fixflat sel |
|---|---|---|
| all | 0.528 / 0.544 | 0.564 / 0.576 |
| EFR_internal | 0.869 / 0.873 | 0.948 / 0.976 |
| external | 0.426 / 0.446 | 0.449 / 0.456 |

The protocol choice does not change any conclusion.

## 3. Results — pooled objects (separate from crop-macro), IoU 0.50

| scope | method | route | filter | GT | pred | TP | FP | FN | P | R | F1 | F1@0.75 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all | frozen | sel | off | 8614 | 6030 | 5123 | 907 | 3491 | 0.850 | 0.595 | 0.700 | 0.663 |
| all | fixflat | sel | off | 8614 | 5643 | 5360 | 283 | 3254 | 0.950 | 0.622 | 0.752 | 0.724 |
| all | frozen | sel | 500 | 8614 | 5571 | 5123 | 448 | 3491 | 0.920 | 0.595 | 0.722 | 0.684 |
| all | fixflat | sel | 500 | 8614 | 5487 | 5360 | 127 | 3254 | 0.977 | 0.622 | 0.760 | 0.732 |
| EFR_internal | frozen | sel | off | 2476 | 2390 | 2261 | 129 | 215 | 0.946 | 0.913 | 0.929 | 0.910 |
| EFR_internal | fixflat | sel | off | 2476 | 2537 | 2423 | 114 | 53 | 0.955 | 0.979 | **0.967** | 0.940 |
| EFR_internal | fixflat | sel | 500 | 2476 | 2447 | 2423 | 24 | 53 | 0.990 | 0.979 | **0.984** | 0.958 |
| external | frozen | sel | off | 6138 | 3640 | 2862 | 778 | 3276 | 0.786 | 0.466 | 0.585 | 0.540 |
| external | fixflat | sel | off | 6138 | 3106 | 2937 | 169 | 3201 | 0.946 | 0.478 | 0.635 | 0.606 |
| external | fixflat | sel | 500 | 6138 | 3040 | 2937 | 103 | 3201 | 0.966 | 0.478 | 0.640 | 0.611 |

Per-specimen pooled rows, ToS route and IoU 0.25/0.75 P/R: `pooled.csv`, `tables.md`.

## 4. Gate activations (E23c: median ≤ 3 ∧ Q75 ≤ 4)

| specimen | frozen: fires / off / crash | fixflat: fires / off / crash |
|---|---|---|
| EFR | 4 / 8 / 0 | **0 / 12 / 0** |
| EFS | 1 / 7 / 0 | **0 / 8 / 0** |
| EMR | 6 / 2 / 0 | 6 / 2 / 0 |
| EMS | 5 / 3 / 0 | 2 / 5 / 1 |
| EOR | 8 / 0 / 0 | 4 / 3 / 1 |
| EOS | 6 / 1 / 1 | 1 / 4 / 3 |

On EFR/EFS the frozen gate fired on 5 boundary crops (specimen fraction 0.41–0.57), the dev crop05 signature; fixflat
removes all 5 (mv 7–47 after the fix) and beats the frozen fallback on 4 of them. Where the frozen gate was off, fixflat is
bit-for-bit the same F1 on all 15 EF crops. Per-crop mv/qv: `per_crop_primary.csv`, `tables.md`.

## 5. Scale: e13 r_est vs granulometric r vs GT inscribed radius (medians per specimen)

| specimen | e13 r_est (min–max) | openrec r | granulometric r_ps_peak | GT inscribed r (diagnostic) |
|---|---|---|---|---|
| EFR | 10.12 (8.0–10.3) | 12–15 | 10.5 | 10.63 |
| EFS | 9.62 (9.0–10.3) | 14–15 | 10.0 | 10.63 |
| EMR | 7.55 (7.0–10.0) | 10–15 | 14.5 | 14.51 |
| EMS | 7.97 (7.35–10.0) | 11–15 | 14.5 | 14.70 |
| EOR | 7.00 (7.0–7.14) | 10–11 | 15.0 | 15.30 |
| EOS | 7.07 (7.0–7.28) | 10–11 | 15.0 | 15.26 |

**e13 underestimates the grain radius by ≈ 2× on the M and O specimens** (ratio 0.46–0.55), while the granulometric
estimate (`batos_v2.scale_granulo`, GT-free, logged only, not used by either method) tracks GT within 0.5 voxel on all
six. The e13 rule is the median of the 1,500 highest EDT local maxima (5³ max filter) of the Otsu foreground; the 1,500 cap
and window are absolute, so on large grains the sample fills with lower medial-ridge maxima. The values sit on
EDT lattice distances (7, √50, 10), i.e. the median lands on ridge voxels, not grain centres.

## 6. Runtime and peak memory (per crop, median [min–max])

| specimen | frozen s | frozen peak MB | fixflat s (end-to-end) | fixflat peak MB |
|---|---|---|---|---|
| EFR | 338 [248–368] | 3312 [3311–3648] | 323 [267–362] | 3164 [3132–3549] |
| EFS | 325 [293–354] | 3312 [3193–3352] | 318 [291–340] | 3204 [3134–3352] |
| EMR | 249 [213–338] | 3501 [3312–3644] | 266 [249–327] | 3364 [3216–3506] |
| EMS | 273 [220–335] | 3504 [3312–3653] | 276 [244–324] | 3355 [3132–3565] |
| EOR | 247 [209–270] | 3565 [3339–3649] | 251 [232–281] | 3453 [3353–3588] |
| EOS | 251 [227–273] | 3511 [3349–3651] | 256 [0–281] | 3370 [3350–3556] |

CPU (AMD Ryzen 7 5700X), runs strictly sequential, another user sweep possibly running concurrently (times are upper
bounds). Fixflat end-to-end = frozen scale + openrec time + fixflat stages. Total wall time 16:03 → 23:42 (incl. a
~10 min stop at the first crash) for 104 runs + evaluation.

## 7. Flags (fixflat F1@0.50 < 0.90) and first-stage forensics

25 of 52 crops are flagged (`flags.csv`). Stage quantities for all 104 runs: `forensics_stages.csv`.

### 7a. The external collapse (21 crops, F1 ≤ 0.21): diverges at stage 03, background rejection

Reference good crop of the same specimen: **EMRGP01_01_test07** (interior, F1 0.941). Failing crop: **EMRGP01_00_test01**
(interior, F1 0.000).

| stage | EMR 01_test07 (good) | EMR 00_test01 (fail) | EFR 00_test01 (dev-like) |
|---|---|---|---|
| 0 e13 r_est → openrec r | 10.0 → 15 | 7.84 → 12 | 10.25 → 15 |
| 1 flatten bg_ref | 201 | 124 | 191 |
| 2 ToS seed components (median area) | 553 (7) | 833 (8) | 602 (15) |
| **3 dark threshold = min(Otsu−10, p40, 135)** | Otsu 114 → **104 (Otsu branch)** | Otsu 145, p40 124 → **124 (p40 branch)** | Otsu 117 → 107 (Otsu) |
| 3 dark fraction of crop / share connected to a face | 0.12 / 0.89 | **0.42 / 0.99** | 0.14 / 0.92 |
| 3 seed components kept | 205 | **0** (all "external_background") | 206 |
| 4 MIN candidates inside internal_dark | 807 | 2 | 483 |
| 5 ToS markers | 258 | 2 | 391 |
| 6 gate mv / qv | 6 / 10 (off) | 2.5 / 2.75 (fires; EDT-Otsu fallback has the same 2 markers) | 45 / 50 |
| F1@0.50 sel | 0.941 | 0.000 | 0.981 |

**Why.** Stage 03 (`src/batos/batos_background.py`, frozen) thresholds the flattened selector at
`min(Otsu − 10, p40, 135)` and rejects every seed that overlaps dark voxels connected to a crop face. On every interior
crop with openrec r ≤ 12, Otsu of the selector rises to 141–147 and the threshold falls to the **40th-percentile branch**,
so by construction ≈ 40 % (up to 83 %) of the crop is "dark", it percolates, ≥ 98 % of it touches a face, and **every
seed is rejected as external background**. Downstream there are 0–13 markers; e11 crashes when there are none.

The branch separates the results completely (fixflat): **p40 branch: 20 crops, F1 0.000–0.209 (mean 0.020); Otsu
branch: 31 crops, F1 0.651–1.000 (mean 0.950).** For interior crops the branch is fully determined by the openrec
radius: r = 15 → Otsu branch 12/12; r ≤ 12 → p40 branch 13/13. Boundary crops (exterior air shifts the histogram) mostly
stay on the Otsu branch under fixflat (19 of 26, r 11–15); the 9 such M/O boundary crops score 0.90–1.00, the 7 M/O
boundary crops on the p40 branch score 0.00–0.21, while all M/O interior crops with r ≤ 12 score 0.

**The same stage also explains the fixflat gain.** Branch per crop, frozen → fixflat: interior crops never change branch
(12 Otsu → Otsu, 13 p40 → p40); on boundary crops the frozen flatten puts **21 of 26 on the p40 branch** and fixflat moves
**14 of them back to the Otsu branch** (5 Otsu → Otsu, 7 p40 → p40). The crop05 support collapse is therefore the same
stage-03 percolation, triggered on boundary crops by the exterior-air bg_ref artefact that fixflat removes, and on M/O
interior crops by the undersized openrec radius, which fixflat does not touch.

**Root-cause chain (observational; no counterfactual run, because no parameter change was allowed):** e13 radius ≈ ½ of
the true radius on M/O grains → openrec r 10–12 instead of ≈ 22 → selector histogram with Otsu ≈ 145 → stage 03 switches to
the p40 branch → percolating "external" dark → all seeds rejected → no markers. The one interior M crop where e13
returned 10.0 (EMR 01_test07, openrec 15) works (0.941); EMS 01_test07 also got openrec 15 and stays on the Otsu branch
(0.651, see 7c). The fixflat change is not involved on interior crops: frozen and fixflat take the same branch on
every one of them.

### 7b. EFR boundary crops (fixflat < 0.90, gate off, Otsu branch)

* **EFRGP01_00_test06** (f 0.41): F1 0.835 = TP 96 / FP 38 / FN 0: every GT grain found; all 38 FP are < 500 voxels
  (filter 500 → F1 **1.000**). Same signature as dev 01_crop05 (small fragments where the specimen meets air).
* **EFRGP01_00_test02** (f 0.56): F1 0.893 = TP 133 / FP 17 / FN 15; filter 500 → 0.924; F1@0.75 0.517 (boundary grains
  segmented with low overlap). Stays on the Otsu branch (threshold 97) but with a low face-connected share (0.73) and
  205 markers for 148 valid GT grains (interior EFR ≈ 395 for ≈ 265); the residual errors are boundary-grain
  under/over-segmentation, not a support collapse (mv 7, qv 13).

### 7c. Interior crops with a low ToS seed count (stage 02 diverges)

* **EFSGP01_01_test05** (0.874, frozen = fixflat) and **EMSGP01_01_test07** (0.651, frozen = fixflat): stage 02 yields only
  171 / 111 seed components with median area 128 / 256 voxels (vs ≈ 600 components of median area ≈ 15 on good interior
  crops); stage 03 keeps nearly all of them, but they are merged cores, so 355 / 186 markers under-seed the crop
  (FN 45 / 74) and the gate statistics drop (mv 12 / 5). Not caused by the flatten fix (identical under frozen).
  Mechanism beyond stage 02 not investigated.

## 8. Conclusions

1. **On the dev material (EF specimens, 20 new crops) the fix generalises.** Fixflat: EFR 0.956 ± 0.048 (filter 500:
   0.985 ± 0.021), EFS 0.968 ± 0.043; gate never fires; pooled F1 0.967 (EFR, 0.984 with filter). Frozen sel reaches only
   0.879 / 0.882 because the crop05 support collapse reappears on 5/10 EF boundary crops. Where frozen did not collapse,
   fixflat is identical.
2. **BA-TOS (frozen and fixflat) does not transfer to the M and O specimens** (crop-macro F1@0.50 sel 0.24–0.48; interior
   crops ≈ 0). The failure is upstream of the flatten fix: the frozen e13 scale estimator underestimates the radius ≈ 2×,
   and the frozen background stage's p40 branch then rejects every seed. Fixflat helps only on those boundary crops that
   escape the p40 branch.
3. Nothing was tuned. Any remedy (e.g. a scale estimator that tracks GT like the granulometric one does, or a
   stage-03 threshold that cannot percolate) is a method change and must be designed on dev data and re-frozen before
   this test set is used again; at that point this test set is no longer untouched for those changes and a further held-out
   set (e.g. other timesteps) should be used.

## Files

`crop_macro.csv`, `pooled.csv`, `per_crop_all_modes.csv` (every crop × method × route × filter × mode × IoU),
`per_crop_primary.csv`, `crop_info.csv` (scale, gate, runtime), `flags.csv`, `forensics_stages.csv`, `tables.md`.
Raw runs: `C:\BA_TOS_DATA\runs\test_protocol\{frozen,fixflat,eval,logs}` (tifs pruned except segmentations/markers).
