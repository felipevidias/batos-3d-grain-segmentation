# Phase 2 — Cellpose baselines vs BA-TOS v3 (dev 12 + held-out test 52 crops)

Status 2026-11-10 (decisive checks 1–6 complete; see the table below). Complete for **cpsam_v2** (main comparison, all analyses) and, for accuracy and runtime only, for
the other three Cellpose models **cpsam (v1), cpdino (DINOv3 ViT-L) and cpdino-vitb (DINOv3 ViT-B)**, all 64 crops,
no failures (§3b). All four Cellpose models agree to within 0.002 F1. v3 = `batos-paper-v3-scalefree` (frozen, tag), BA-TOS_sel route; the gate never fires on
these crops, so sel = ToS. Nothing was tuned on GT. No size filter unless stated.

Sources of record: `tables.md` (crop-macro and pooled, filter none/500, runtime), `protocol_check/tables.md`
(no-marker vs strict), `shape/tables.md` (morphology, boundary, scoring modes, contacts),
`robustness/tables.md`, `synthetic_scale/tables.md`. Rows `cellpose_<tag>/<crop>` in `results/all_runs.csv`.

**Decisive checks (request of 2026-11-10) → where answered:**

| # | check | section | one-line answer |
|---|---|---|---|
| 1 | leakage | §2 | No micro-CT / granular / sand / soil / Zenodo 8014905 data in the documented Cellpose-SAM set (confirmed in the preprint's text). The June 2026 models (cpsam_v2, cpdino*) have no published dataset list. |
| 2 | strict mode, 64 crops | §4 | Strict reverses the ranking (test 0.586 vs 0.776), only because cpsam_v2 also segments the face-truncated grains that are excluded from GT; counting them as GT gives 0.979 vs 0.604. |
| 3 | per-grain morphology | §5 | cpsam_v2 is closer to GT shape (sphericity, elongation, centroid, MSD, HD95) on every specimen; its one bias is a ~1-voxel dilation. |
| 4 | contacts | §6 | Contact recall vs GT 0.998 (cpsam_v2) vs 0.801 (v3); vs the authors' ContactPoint lists 1.000 vs 0.782. |
| 5 | robustness (crop04, EORGP01_00_test01, + test04) | §7 | v3 collapses at 5 % noise on all three crops; cpsam_v2 ≥ 0.977 everywhere. |
| 6 | synthetic scale 0.7/1.4/2.0 | §7b | cpsam_v2 1.000 / 0.999 / 0.963 (1.000 at every factor with the size filter); v3 0.972 / 0.979 / 0.900. |

## 1. Setup

| item | value |
|---|---|
| env | `C:\BA_TOS_DATA\envs\cellpose`, Python 3.12.10, torch 2.11.0+cu128, cellpose 4.2.1.1, dinov3 @ 6876159a (no-deps); `environments/cellpose.txt` |
| GPU | RTX 5060 8 GB (sm_120), CUDA available, conv3d GPU = CPU within TF32 (9e-4) |
| models offered (`cellpose.models.MODEL_NAMES`) | cpsam_v2 (default, June 2026), cpdino (DINOv3-ViTL, June 2026), cpdino-vitb (DINOv3-ViTB), cpsam (original, April 2025) |
| config (one per model, no per-crop tuning) | `configs/cellpose_<model>.yaml`: do_3D, z_axis 0, diameter None (cellpose 4 has no size model: no rescaling), flow3D_smooth 0, all else cellpose defaults (resolved signature in every run.json) |
| runner / evaluator | `scripts/run_cellpose.py` (prediction only, GT never read); `scripts/evaluate_cellpose.py` (evaluate_v2, no-marker fraction 0.50, IoU .25/.50/.75, filter none/500; v3 rows read from the frozen v3 evaluations in the same mode) |

## 2. Leakage check (training data of Cellpose-SAM / CellposeDINO)

- **cpsam (original Cellpose-SAM, April 2025). Confirmed from the preprint's full text** (bioRxiv
  10.1101/2025.04.28.651001v1, posted 1 May 2025, PDF read 2026-11-10). Verbatim: *"The model was then trained on an
  updated dataset of cells and nuclei containing 22,826 train images with a combined 3,341,254 training ROIs. This
  dataset combines major currently available datasets: Cellpose, Cellpose Nuclei, Omnipose, TissueNet, LiveCell,
  YeaZ, DeepBacs, Neurips 2022, MoNuSeg, MoNuSAC, CryoNuSeg, NuInsSeg, BCCD, CPM 15+17, TNBC, LynSec, IHC TMA,
  CoNIC, PanNuke"*. The encoder is initialised from SAM, pretrained on SA-1B (natural images). Training is **2D only**.
  3D masks come from averaging 2D flows over XY/XZ/YZ slices at test time. The only 3D data in the paper
  (BlastoSPIM, PlantSeg) is used for finetuning experiments and was explicitly *not* used for training. The data
  availability statement points to cellpose.org/dataset and to the cited labs' datasets. None of the datasets is
  micro-CT, X-ray, granular, sand or soil, and Zenodo 8014905 / Pinzón is not listed.
- **Relevant to §7 (not leakage).** The preprint's training augmentations add, on 50 % of every batch, *"Poisson
  noise, Gaussian blurring, downsampling, and anisotropic blurring with downsampling"*, plus random rescaling by
  0.25–4× around a 30-px mean diameter. Robustness to noise, downsampling and scale was trained in explicitly.
  The §7 and §7b results are therefore expected behaviour, not luck.
- **The one non-cell component** is inherited from the original Cellpose dataset (cyto/cyto2; Stringer et al., Nat.
  Methods 2021). It includes a small set (n = 98) of *non-microscopy* images "which contained large numbers of
  repeated objects such as fruits, rocks and jellyfish", collected by internet search. These are 2D natural photographs,
  not CT. Packed rocks are the visually closest content to our data. This is not label leakage of our data, but it is
  prior exposure to the "packed convex objects" appearance.
- **cpsam_v2, cpdino, cpdino-vitb (June 2026).** The Cellpose model docs give only the backbone, the release date, a
  "fix in the training for low contrast regions" (cpsam_v2), and the training size range (diameters 7.5–120 px, mean
  30). **No dataset list was found** for these three; I could not verify their training data. The DINOv3 backbone is
  self-supervised pre-trained by Meta on LVD-1689M (web images, no segmentation labels). It may contain photos of sand
  or gravel, but not our labels.
- **cpsam_v2 / cpdino / cpdino-vitb, re-checked 2026-11-10.** The GitHub README (MouseLand/cellpose) says only that
  cpsam_v2 *"predicts fewer spurious masks in low-contrast regions"* and that cpdino models are *"based on DINOv3"*. It
  cites no paper for these models. *"All Cellpose models are trained on data that is licensed under CC-BY-NC."*
  The docs page (models.html) gives the backbone, the date and the 7.5–120 px diameter range only. No model card or
  dataset list exists; to be asked of the authors before submission.
- **Verdict.** There is no evidence of micro-CT, granular-material or Zenodo 8014905 data in any documented Cellpose
  training set; the Cellpose-SAM list is now confirmed from the primary text. Residual uncertainty is limited to the
  undocumented June 2026 training sets. The original cpsam (whose list is documented) gives the same results
  (§3b: test 0.999), so even a hypothetical addition of our data to the 2026 sets could not explain the result.
- **A separate open point, not leakage.** How the Pinzón GT labels were produced (an automatic pipeline such as
  thresholding + watershed in spam, or manual) is not documented in this repo. If the GT is itself a watershed of a
  clean binarisation, "agreement with GT" partly measures agreement with that pipeline.

Sources: bioRxiv 10.1101/2025.04.28.651001v1 full text (Cellpose-SAM, methods "Model design" and "Whole-dataset
training"); github.com/MouseLand/cellpose README; cellpose.readthedocs.io/en/latest/models.html; Stringer
et al. 2021 Nat. Methods / bioRxiv 10.1101/2020.02.02.931238 (non-microscopy subset); Janelia figshare 13270466
(Cellpose training dataset); ai.meta.com/dinov3 (DINOv3).

## 3. Main comparison (no-marker mode, no filter)

Crop-macro F1@0.50 (mean ± SD, min) and pooled F1@0.50:

| scope | n | cpsam_v2 macro | v3 macro | cpsam_v2 pooled (TP/FP/FN) | v3 pooled (TP/FP/FN) |
|---|---|---|---|---|---|
| dev (EFR) | 12 | **0.998** ± 0.004 (0.987) | 0.979 ± 0.011 (0.958) | 0.998 (2808/13/0) | 0.978 (2711/25/97) |
| test all | 52 | **0.998** ± 0.004 (0.980) | 0.896 ± 0.168 (0.398) | 0.998 (8613/32/1) | 0.907 (7437/354/1177) |
| test EFR | 12 | 0.999 | 0.977 | 0.999 | 0.976 |
| test EFS | 8 | 0.999 | 0.964 | 0.999 | 0.960 |
| test EMR | 8 | 0.998 | 0.966 | 0.998 | 0.959 |
| test EMS | 8 | 0.999 | 0.970 | 0.998 | 0.965 |
| test EOR | 8 | 0.996 | 0.735 | 0.996 | 0.714 |
| test EOS | 8 | 0.998 | 0.727 | 0.998 | 0.709 |

- F1@0.75 test: cpsam_v2 0.997, v3 0.853 (macro). Filter 500: cpsam_v2 dev 1.000; full tables in `tables.md`.
- No cpsam_v2 crop is below 0.90 (v3: 10 of 52, all EO interior plus EFSGP01_01_test05 and EOSGP01_01_test08).
- An independent re-implementation of the IoU on v3's worst crop (EORGP01_01_test05) confirms 162/162 valid GT
  matched by cpsam_v2 at IoU ≥ 0.75 (median 0.913). The prediction is not a copy of the GT (foreground 0.711 vs 0.661).

## 3b. All four Cellpose models (accuracy, same fixed config, no filter)

| model | dev macro F1@0.50 (min) | test macro F1@0.50 ± SD (min) | test pooled TP/FP/FN | test pooled F1@0.50 / @0.75 | test mIoU@0.50 | EO macro (EOR / EOS) |
|---|---|---|---|---|---|---|
| cpsam_v2 | 0.998 (0.987) | 0.998 ± 0.004 (0.980) | 8613 / 32 / 1 | 0.998 / 0.997 | 0.899 | 0.996 / 0.998 |
| cpsam (v1) | 1.000 (0.998) | 0.999 ± 0.002 (0.988) | 8613 / 17 / 1 | 0.999 / 0.998 | 0.912 | 0.998 / 1.000 |
| cpdino | 0.999 (0.994) | 0.999 ± 0.002 (0.988) | 8613 / 12 / 1 | 0.999 / 0.999 | **0.928** | 0.998 / 1.000 |
| cpdino-vitb | 0.999 (0.995) | 0.998 ± 0.004 (0.977) | 8613 / 18 / 1 | 0.999 / 0.998 | 0.910 | 0.995 / 1.000 |
| v3 | 0.979 (0.958) | 0.896 ± 0.168 (0.398) | 7437 / 354 / 1177 | 0.907 / 0.875 | 0.902 | 0.735 / 0.727 |

- Detection is saturated for every model: all four miss the same single GT grain on test (FN 1) and differ only by
  12–32 FPs over 8,614 grains. The models differ only in boundary overlap: cpdino has the highest matched IoU
  (0.928 test, 0.912 dev), above v3 (0.902 / 0.916). cpsam_v2, the default, has the lowest (0.899 / 0.888).
  This is consistent with the ~1-voxel dilation found for cpsam_v2 in §5. Shape/contact analyses (§5–6) were run
  for cpsam_v2 only, so they describe the *weakest-boundary* Cellpose model.
- The conclusion is about the model family, not about one checkpoint.

## 3c. Sensitivity of cpsam_v2 (dev 12 only, one fixed change at a time, no filter)

| run | macro F1@0.50 ± SD (min) | pooled TP/FP/FN | F1@0.75 macro | mIoU@0.50 | s/crop |
|---|---|---|---|---|---|
| baseline (flow3D_smooth 0, diameter None) | 0.998 ± 0.004 (0.987) | 2808 / 13 / 0 | 0.998 | 0.888 | 30 |
| flow3D_smooth 1 | 0.998 ± 0.003 (0.989) | 2808 / 10 / 0 | 0.998 | 0.888 | 30 |
| flow3D_smooth 2 | 0.999 ± 0.002 (0.994) | 2808 / 7 / 0 | 0.999 | 0.888 | 31 |
| diameter = 2·r_granulo (GT-free, = 20 on every dev crop → rescale 1.5) | 0.997 ± 0.004 (0.987) | 2808 / 19 / 0 | 0.997 | 0.861 | 171 |

- Every setting keeps recall at 1.000. Only the number of small spurious FPs (7–19) changes, and the size filter 500
  removes all of them (1.000 for every run). The result is insensitive to the two documented knobs, so the default
  config is not a cherry-picked operating point.
- flow3D_smooth makes no difference to boundaries (mIoU 0.888 in every run). Upscaling with the BA-TOS scale estimate
  makes boundaries worse (mIoU 0.861; cause not measured) and is 5.7× slower. Cellpose's own
  scale (no rescale) is the right choice here. Test was not run with these settings: nothing was selected.
- Tables: `sensitivity/tables.md`; rows `cellpose_cpsam_v2__<override>/<crop>` in `all_runs.csv`.

## 4. Protocol check — what the ignore rule does

crop04: cpsam_v2 630 predictions, **358 ignored** by the no-marker rule, 272 evaluable, 0 FP, pred/GT 1.000. v3 396
predictions, 131 ignored, 265 evaluable, 2 FP, pred/GT 0.974. The crop has 377 excluded GT grains (face-touching or
< 1000 voxels).

Four scoring modes, pooled F1@0.50 (TP/FP/FN):

| scope | method | no-marker (paper) | strict: nothing ignored | all GT incl. face-touching ≥ 1000 | all GT, every label |
|---|---|---|---|---|---|
| dev | cpsam_v2 | 0.998 (2808/13/0) | 0.615 (2808/3519/0) | 0.951 (5739/588/0) | **0.979** (6308/23/243) |
| dev | v3 | 0.978 (2711/25/97) | **0.829** (2773/1107/35) | 0.744 (3654/424/2085) | 0.687 (3650/428/2901) |
| test 52 | cpsam_v2 | 0.998 (8613/32/1) | 0.586 (8613/12193/1) | 0.952 (18900/1908/3) | **0.979** (20740/98/780) |
| test 52 | v3 | 0.907 (7437/354/1177) | **0.776** (7697/3516/917) | 0.657 (10064/1691/8839) | 0.604 (10056/1699/11464) |

v3 under its own **seed-ignore** rule (paper §III-C, the SIBGRAPI protocol; not applicable to Cellpose, which has
no seeds): dev pooled 0.988 (2765/25/43), crop-macro 0.989; test pooled 0.921 (7652/358/962), crop-macro 0.912. So
on all 64 crops (pooled F1@0.50, dev / test): v3 seed-ignore 0.988 / 0.921, no-marker 0.978 / 0.907, strict
0.829 / 0.776; cpsam_v2 no-marker 0.998 / 0.998, strict 0.615 / 0.586. Even v3's most favourable protocol is below
cpsam_v2's paper-protocol score on both splits.

Interpretation: strict mode (nothing ignored, standard exclusions) reverses the ranking. But its extra "FPs" are
predictions of the **excluded** GT grains: the grains truncated by the crop faces (and small fragments). cpsam_v2
segments them; v3 rejects boundary-connected candidates by design. When those grains are counted as GT ("all GT,
every label"), cpsam_v2 recovers them almost completely (20,740 TP, 98 FP, F1 0.979). v3 does not (0.604). Strict
mode therefore measures how much a method refrains from segmenting excluded grains, not segmentation quality. The
paper protocol (exclusion + ignore rule) and the all-GT mode agree on the ranking. Per-specimen strict numbers:
`protocol_check/tables.md`.

## 5. Shape and boundary fidelity (matched at IoU 0.50, paired on GT grains matched by both)

Mean signed error [median |error|], voxel units, test 52 (dev in `shape/tables.md`):

| scope | method | n | volume % | d_eq | surface % | sphericity | elongation | centroid | MSD | HD95 |
|---|---|---|---|---|---|---|---|---|---|---|
| test all | cpsam_v2 | 7436 | +9.8 [9.9] | +0.96 [0.96] | +5.2 [5.3] | +0.010 [0.009] | −0.012 [0.015] | 0.18 [0.15] | **0.64** [0.64] | **1.08** [1.00] |
| test all | v3 | 7436 | +12.8 [6.4] | +1.12 [0.63] | +27.2 [7.6] | −0.106 [0.048] | +0.175 [0.017] | 2.13 [0.11] | 1.43 [0.49] | 6.51 [1.41] |
| test EFR | cpsam_v2 | 2394 | +11.0 [10.9] | +1.01 [1.01] | +6.0 [6.0] | +0.010 [0.009] | −0.025 [0.025] | 0.18 [0.15] | 0.64 [0.64] | 1.06 [1.00] |
| test EFR | v3 | 2394 | +5.8 [6.5] | +0.44 [0.64] | +17.7 [2.8] | −0.082 [0.030] | +0.078 [0.015] | 1.50 [0.08] | 1.09 [0.48] | 4.83 [1.00] |
| test EO (R) | cpsam_v2 | 691 | +7.9 [7.9] | +0.84 [0.84] | +3.9 [3.9] | +0.011 [0.011] | +0.005 [0.008] | 0.16 [0.13] | 0.59 [0.59] | 1.06 [1.00] |
| test EO (R) | v3 | 691 | +28.3 [15.7] | +2.59 [1.62] | +48.0 [36.0] | −0.158 [0.145] | +0.341 [0.048] | 3.34 [2.16] | 2.11 [1.31] | 9.72 [5.75] |

Mean descriptors on the paired test set: GT sphericity 0.882 → cpsam_v2 0.892, v3 0.777; GT elongation 1.348 →
cpsam_v2 1.336, v3 1.523.

**Which is closer to GT shape?**
- **cpsam_v2** is closer on sphericity, elongation, centroid, MSD and HD95, on every specimen. Its error is one
  systematic bias: a ≈ 1-voxel outward boundary (volume +8 to +11 %, d_eq ≈ +1 voxel, HD95 ≈ 1.0, very low spread).
- **v3** has the smaller *median* volume and d_eq error on the EF specimens (no dilation bias). Its mean errors are
  dominated by a tail of badly shaped grains (mean HD95 4.6–10 voxels; mean centroid error 1.4–3.6 voxels vs
  median 0.08–2), and its sphericity is biased low (−0.08 to −0.17): rough watershed boundaries and mis-split grains.
- The dilation bias of cpsam_v2 looks correctable (e.g. a 1-voxel erosion). That would be a post-hoc
  calibration, so it was **not** applied. Testing it would have to be designed on dev and frozen.

## 6. Contacts (6-connectivity adjacency between two valid GT grains)

| scope | method | vs GT 6-adj: recall / precision | vs authors' ContactPoint: recall / precision |
|---|---|---|---|
| dev | GT 6-adj itself | — | 1.000 / 0.618 |
| dev | cpsam_v2 | 0.999 / 0.962 | 1.000 / 0.595 |
| dev | v3 | 0.938 / 0.967 | 0.951 / 0.606 |
| test 52 | GT 6-adj itself | — | 1.000 / 0.668 |
| test 52 | cpsam_v2 | 0.998 / 0.962 | 1.000 / 0.644 |
| test 52 | v3 | 0.801 / 0.930 | 0.782 / 0.607 |
| test EO (R / S) | cpsam_v2 | 0.994 / 0.970, 0.998 / 0.964 | 0.996 / 0.667, 1.000 / 0.679 |
| test EO (R / S) | v3 | 0.467 / 0.829, 0.440 / 0.813 | 0.396 / 0.483, 0.382 / 0.495 |

- The authors' contact lists (Contacts/ContactPoint/<scan>.txt; available for all six specimens, 64/64 crops) use
  full-volume label ids that match the crop GT: every listed contact point between two valid grains of a crop lies
  inside the crop.
- Every authors' contact is a 6-adjacent GT pair (recall 1.000). However, 6-adjacency finds ~1.5× more pairs
  (precision 0.62–0.67), so the authors apply a stricter contact criterion (probably a minimum contact area or a
  different detection). The precision of both methods against the authors' list is therefore capped near 0.62–0.67.
  Recall against the authors' list is the meaningful number there.
- cpsam_v2 recovers essentially all contacts. v3 loses contacts in proportion to its missed and mis-split grains (EO
  recall ≈ 0.4).

## 7. Robustness (degraded inputs; frozen methods, no retuning)

`scripts/robustness_cellpose_v3.py`. Crops: EFRGP01_00_crop04 (dev), EORGP01_00_test01 (requested; an EO interior
crop on which v3 already fails when clean, 0.544, and cpsam_v2 has its lowest test score, 0.980) and
EORGP01_00_test04 (v3's best EO crop, clean 0.985, so degradation is visible).
Additive Gaussian noise, sigma = 2 / 5 / 10 % of the grain-air contrast (difference of the Otsu class means,
GT-free; seed 20261109). 0.5× = 2×2×2 block-mean downsampling to 100³; predictions upsampled ×2 (nearest) and
scored on the original GT. Both methods frozen. v3 adapts its scale itself (r_granulo 10 → 5 and 15 → 7 at 0.5×).

No-marker mode, F1@0.50 (F1@0.75):

| crop | method | clean | noise 2 % | noise 5 % | noise 10 % | 0.5× |
|---|---|---|---|---|---|---|
| EFRGP01_00_crop04 | cpsam_v2 | 1.000 (1.000) | 1.000 (1.000) | 1.000 (1.000) | 1.000 (1.000) | 1.000 (1.000) |
| EFRGP01_00_crop04 | v3 | 0.980 (0.972) | 0.912 (0.908) | 0.127 (0.103) | 0.000 (0.000) | 0.662 (0.662) |
| EORGP01_00_test01 | cpsam_v2 | 0.980 (0.980) | 0.980 (0.980) | 0.977 (0.977) | 0.977 (0.977) | 0.985 (0.985) |
| EORGP01_00_test01 | v3 | 0.544 (0.471) | 0.685 (0.555) | 0.062 (0.052) | 0.000 (0.000) | 0.979 (0.967) |
| EORGP01_00_test04 | cpsam_v2 | 1.000 (1.000) | 1.000 (1.000) | 1.000 (1.000) | 0.990 (0.990) | 1.000 (1.000) |
| EORGP01_00_test04 | v3 | 0.985 (0.917) | 0.941 (0.861) | 0.000 (0.000) | 0.000 (0.000) | 0.976 (0.812) |

- The clean v3 re-runs reproduce the frozen scores (0.980, 0.544, 0.985), so the robustness harness is faithful.
- **EORGP01_00_test01** confirms the capacity mechanism from both sides. Clean, v3 keeps only 156 ToS markers for
  170 valid grains (96 FN). At 10 % noise no ToS marker survives (EMPTY, scored as an empty segmentation). At 0.5×
  downsampling the grains shrink back under the candidate cap and v3 jumps to **0.979** (234 markers). cpsam_v2 is
  flat (0.977–0.985); its 6–7 FPs and 1 FN are the same in every condition.
- **Native noise level.** The intrinsic noise of the scans is about 4.5 % (crop04) and 1.9 % (EO crop) of contrast.
  This is a robust MAD estimate of neighbour differences and slightly overestimated because it includes edges.
  Added noise of 5 % therefore roughly doubles the noise variance on crop04.
- **v3 failure mechanism under noise.** Noise multiplies small dark extrema, which saturate the ToS candidate
  capacity. The surviving ToS markers fall from 396 → 374 → 78 → 5 (crop04) and from 149 → 143 → 3 → none (EO).
  The support-collapse gate fires only at 10 % on crop04 (mv 2), and its EDT-Otsu fallback does not rescue the
  result (5 predictions). EORGP01_00_test04 at 10 % gives no segmentation at all. It is the same capacity limit as
  the EO failure in `results/test_v3/REPORT.md`, now triggered by noise instead of grain size.
- **Downsampling.** v3 over-segments crop04 at 0.5× (691 markers for ~270 valid grains; F1 0.662) but holds on EO (0.976).
- **cpsam_v2** is unaffected except 0.990 at 10 % noise on the EO crop. Strict-mode numbers are in
  `robustness/tables.md`. They show the same pattern; cpsam_v2's strict ~0.60 is the face-grain effect of §4.
- Strict-mode collapse of v3 (0.149 / 0.000 at 5 %) confirms that the noise failure is not an artefact of the ignore rule.

## 7b. Synthetic scale set (factors 0.7 / 1.4 / 2.0 of the 12 dev crops)

`scripts/evaluate_cellpose_synth.py`; crops from `scripts/make_synthetic_scale.py` (dev crops resampled about their
centre; valid GT ≥ round(1000 f³) voxels). cpsam_v2 with the same fixed config (no diameter); v3 rows from the v3
scale test. No-marker mode, size filter none; F1@0.50 crop-macro mean ± SD (min) and pooled TP/FP/FN.

| factor | cpsam_v2 macro | cpsam_v2 pooled | v3 macro | v3 pooled |
|---|---|---|---|---|
| 0.7 | **1.000** ± 0.000 (0.998) | 8939 / 5 / 0 | 0.972 ± 0.012 (0.946) | 8813 / 342 / 126 |
| 1.0 (original) | **0.998** ± 0.004 (0.987) | 2808 / 13 / 0 | 0.979 ± 0.011 (0.958) | 2711 / 25 / 97 |
| 1.4 | **0.999** ± 0.002 (0.994) | 877 / 1 / 0 | 0.979 ± 0.007 (0.962) | 849 / 9 / 28 |
| 2.0 | **0.963** ± 0.028 (0.919) | 188 / 14 / 0 | 0.900 ± 0.149 (0.571) | 174 / 21 / 14 |

- With the size filter (500 f³ for cpsam_v2, v3's own 500 s³): cpsam_v2 **1.000 at every factor** (0 FP, 0 FN);
  v3 0.991 / 0.981 / 0.984 / 0.960.
- cpsam_v2's only scale effect is at 2×: a few small spurious fragments (14 FPs over 12 crops), never a missed grain.
  v3 loses grains at 2× (14 FN; worst crop 0.571, EFRGP01_01_crop05_x2).
- mIoU@0.50: cpsam_v2 0.870 / 0.888 / 0.877 / 0.891, v3 0.885 / 0.916 / 0.904 / 0.914. v3's matched grains overlap GT
  slightly better at every scale (no dilation bias). This is the one metric where v3 leads, on the grains it finds.
- The 0.7–2× range stays within Cellpose's trained diameter range (7.5–120 px, with 0.25–4× rescaling augmentation,
  §2). This test therefore cannot find a scale regime outside Cellpose's training.

## 8. Runtime and memory (200³ crop)

| method | device | median s/crop [range] | peak VRAM | peak host RAM |
|---|---|---|---|---|
| cpsam_v2 | RTX 5060 (GPU) | 30 [29–46] | 1.8 GB reserved | 3.5 GB |
| cpsam (v1) | RTX 5060 (GPU) | 31 [30–60] | 1.8 GB reserved | 3.5 GB |
| cpdino (ViT-L) | RTX 5060 (GPU) | 58 [57–77] | 1.7 GB reserved | 3.3 GB |
| cpdino-vitb | RTX 5060 (GPU) | 23 [22–28] | 1.2 GB reserved | 2.5 GB |
| v3 | Ryzen 7 5700X (CPU) | 397 dev / 397 test | — | 3.4 GB |

Cellpose time is `model.eval` only (model load excluded, not measured). v3 is the whole pipeline, CPU, not
contention-controlled. These are different devices; this is not a like-for-like speed comparison.

## 9. Implications for the paper (one page)

**What the data say.** Run out of the box (one fixed config per model, no fine-tuning, no per-crop tuning, no GT
use), every Cellpose 4 model tested (cpsam_v2, cpsam, cpdino, cpdino-vitb) segments these micro-CT grain crops
almost perfectly under the paper's protocol. Test crop-macro F1@0.50 is 0.998–0.999 over 52 held-out crops from six
specimens, against 0.896 for frozen BA-TOS v3 (0.912 under v3's own seed-ignore rule). The gap is not a protocol
artefact (§4: strict mode reverses the ranking only by counting segmented face-truncated grains as FP; counting them
as GT gives 0.979 vs 0.604). It is not leakage (§2: documented training data is 2D cells/nuclei; the original cpsam,
whose list is documented, gives the same result). It is not one checkpoint (§3b) or one setting (§3c). It holds for
shape (§5), contacts (§6), noise and resolution (§7) and scale 0.7–2× (§7b). Where BA-TOS fails (EO specimens,
noise ≥ 5 %), the failure is one mechanism, ToS candidate capacity, which the data confirm from both sides (§7,
EORGP01_00_test01).

**What cannot be claimed.**
- That BA-TOS is more accurate than learned / foundation models on micro-CT grains, on any metric reported here
  except mean matched IoU of the grains it does find (v3 0.902–0.916 vs cpsam_v2 0.888–0.899; cpdino 0.928 beats it).
- That BA-TOS is more robust, or that a "harder regime" (noise, low resolution, other scales) favours it: §7 and §7b
  show the opposite, and Cellpose-SAM was explicitly trained with noise, blur, downsampling and 0.25–4× rescaling.
- Anything based on the strict mode.

**What can be claimed honestly.**
1. A training-free, GPU-free, fully inspectable pipeline: every marker traces to a ToS node with a score, failures
   are diagnosable and GT-free detectable (capacity and support-collapse statistics). Real, but secondary next to a
   0.10 F1 gap on held-out data.
2. Boundary placement without a dilation bias on the grains it finds (lower median volume error on EF; higher matched
   IoU than cpsam_v2). Narrow; a 1-voxel calibration of Cellpose would probably close it, and cpdino already does.
3. Licensing and provenance: no CC-BY-NC weights, no DINOv3 licence. Relevant to some users; not a scientific result.

**Recommended framing (decision for the authors).**
- **(A) Benchmark and protocol paper, recommended.** Contribution: a GT-free-defined 64-crop benchmark over six
  specimens; an analysis of how face exclusion and the ignore rule move rankings (§4, the cleanest new result); shape,
  contact and robustness evaluation against the authors' own contact data; BA-TOS v3 as the transparent morphological
  reference with its failure mechanism explained, Cellpose 4 as the learned reference. Everything needed is done.
- **(B) Hybrid.** Viable only if it beats cpsam_v2/cpdino on boundary, shape or contact metrics; F1@0.50 has no room
  left (≤ 0.002). Candidates: a ToS-based boundary refinement of Cellpose masks, or GT-free QC that flags failures. It
  must be designed on dev and frozen before test, and §5 should first be rerun with cpdino as the stronger baseline.
- **(C) Harder regime: not supported.** Needs data outside Cellpose's training experience (e.g. fractured,
  non-convex or very large grains, other modalities). None is in hand.

**Before submission.** Ask the Cellpose authors for the June 2026 training sets. Document how the Pinzón GT was
produced: if it is itself a watershed of a clean threshold, "agreement with GT" partly measures agreement with that
pipeline, for both methods. Report v3's mIoU advantage and cpsam_v2's dilation bias explicitly, rather than only F1.
