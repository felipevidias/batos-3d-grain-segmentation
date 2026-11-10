# Phase 2 — Cellpose baselines vs BA-TOS v3 (dev 12 + held-out test 52 crops)

Status 2026-11-10. Complete for **cpsam_v2** (main comparison, all analyses) and, for accuracy and runtime only, for
the other three Cellpose models **cpsam (v1), cpdino (DINOv3 ViT-L) and cpdino-vitb (DINOv3 ViT-B)**, all 64 crops,
no failures (§3b). All four Cellpose models agree to within 0.002 F1. v3 = `batos-paper-v3-scalefree` (frozen, tag), BA-TOS_sel route; the gate never fires on
these crops, so sel = ToS. Nothing was tuned on GT. No size filter unless stated.

Sources of record: `tables.md` (crop-macro and pooled, filter none/500, runtime), `protocol_check/tables.md`
(no-marker vs strict), `shape/tables.md` (morphology, boundary, scoring modes, contacts),
`robustness/tables.md`. Rows `cellpose_cpsam_v2/<crop>` in `results/all_runs.csv`.

## 1. Setup

| item | value |
|---|---|
| env | `C:\BA_TOS_DATA\envs\cellpose`, Python 3.12.10, torch 2.11.0+cu128, cellpose 4.2.1.1, dinov3 @ 6876159a (no-deps); `environments/cellpose.txt` |
| GPU | RTX 5060 8 GB (sm_120), CUDA available, conv3d GPU = CPU within TF32 (9e-4) |
| models offered (`cellpose.models.MODEL_NAMES`) | cpsam_v2 (default, June 2026), cpdino (DINOv3-ViTL, June 2026), cpdino-vitb (DINOv3-ViTB), cpsam (original, April 2025) |
| config (one per model, no per-crop tuning) | `configs/cellpose_<model>.yaml`: do_3D, z_axis 0, diameter None (cellpose 4 has no size model: no rescaling), flow3D_smooth 0, all else cellpose defaults (resolved signature in every run.json) |
| runner / evaluator | `scripts/run_cellpose.py` (prediction only, GT never read); `scripts/evaluate_cellpose.py` (evaluate_v2, no-marker fraction 0.50, IoU .25/.50/.75, filter none/500; v3 rows read from the frozen v3 evaluations in the same mode) |

## 2. Leakage check (training data of Cellpose-SAM / CellposeDINO)

- **cpsam (original Cellpose-SAM, April 2025).** The preprint (Pachitariu, Rariden, Stringer, bioRxiv
  10.1101/2025.04.28.651001) trains on 22,826 images / 3,341,254 ROIs of **cells and nuclei**. Named sources: Cellpose
  (cyto2), Cellpose Nuclei, Omnipose, TissueNet, LiveCell, YeaZ, DeepBacs, NeurIPS 2022, MoNuSeg, MoNuSAC, CryoNuSeg,
  NuInsSeg, BCCD, CPM 15+17, TNBC, LynSec, IHC TMA, CoNIC, PanNuke. These are all 2D light-microscopy or histology.
  None is micro-CT, granular, sand or soil, and Zenodo 8014905 / Pinzón is not listed.
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
- **Verdict.** There is no evidence of micro-CT, granular-material or Zenodo 8014905 data in any documented Cellpose
  training set. Residual uncertainty: (a) the June 2026 models' data are undocumented; (b) the bioRxiv full text
  returned HTTP 429, so the dataset list above comes from search excerpts of the preprint. It should be confirmed
  against the preprint's methods section before it is cited.
- **A separate open point, not leakage.** How the Pinzón GT labels were produced (an automatic pipeline such as
  thresholding + watershed in spam, or manual) is not documented in this repo. If the GT is itself a watershed of a
  clean binarisation, "agreement with GT" partly measures agreement with that pipeline.

Sources: bioRxiv 10.1101/2025.04.28.651001 (Cellpose-SAM); cellpose.readthedocs.io/en/latest/models.html; Stringer
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

`scripts/robustness_cellpose_v3.py`. Crops: EFRGP01_00_crop04 (dev) and EORGP01_00_test04. The EO crop is a
boundary crop and v3's best EO crop (clean 0.985); every EO interior crop is already a v3 capacity failure.
Additive Gaussian noise, sigma = 2 / 5 / 10 % of the grain-air contrast (difference of the Otsu class means,
GT-free; seed 20261109). 0.5× = 2×2×2 block-mean downsampling to 100³; predictions upsampled ×2 (nearest) and
scored on the original GT. Both methods frozen. v3 adapts its scale itself (r_granulo 10 → 5 and 15 → 7 at 0.5×).

No-marker mode, F1@0.50 (F1@0.75):

| crop | method | clean | noise 2 % | noise 5 % | noise 10 % | 0.5× |
|---|---|---|---|---|---|---|
| EFRGP01_00_crop04 | cpsam_v2 | 1.000 (1.000) | 1.000 (1.000) | 1.000 (1.000) | 1.000 (1.000) | 1.000 (1.000) |
| EFRGP01_00_crop04 | v3 | 0.980 (0.972) | 0.912 (0.908) | 0.127 (0.103) | 0.000 (0.000) | 0.662 (0.662) |
| EORGP01_00_test04 | cpsam_v2 | 1.000 (1.000) | 1.000 (1.000) | 1.000 (1.000) | 0.990 (0.990) | 1.000 (1.000) |
| EORGP01_00_test04 | v3 | 0.985 (0.917) | 0.941 (0.861) | 0.000 (0.000) | 0.000 (0.000) | 0.976 (0.812) |

- The clean v3 re-runs reproduce the frozen scores (0.980 and 0.985), so the robustness harness is faithful.
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

## 9. Implications for the journal paper

**The accuracy claim is gone.** Used out of the box (default settings, one fixed config, no per-crop tuning, no
fine-tuning), Cellpose-SAM v2 segments these micro-CT grain crops essentially perfectly under the paper's protocol:
0.998 crop-macro F1@0.50 on 52 held-out crops from six specimens, minimum 0.980. It beats BA-TOS v3 on every specimen,
on F1@0.75, on shape descriptors, on boundary distance and on contacts. The gap is largest exactly where v3 is
weakest (EO, 0.73 → 0.997). The two scoring modes that change the treatment of excluded grains do not rescue
BA-TOS: the "strict" mode favours v3 only because it counts segmented face-truncated grains as FPs, and counting
those grains as GT gives cpsam_v2 0.979 vs v3 0.604. A journal paper whose thesis is "BA-TOS outperforms learned
models on micro-CT grains" would be refuted by its own baseline. The robustness test widens the gap: v3 collapses at
5 % added noise (≈ 2× the native noise variance), while cpsam_v2 stays at ≥ 0.990. The original Cellpose-SAM
checkpoint (cpsam) and both CellposeDINO models (cpdino 0.999, cpdino-vitb 0.998) give the same picture, so this
is a property of the model family, not of one lucky checkpoint.

**What BA-TOS still has, honestly stated:**
1. *No learning, no GPU, no external weights, fully inspectable.* Every marker traces to a Tree-of-Shapes node with
   a score; the failure mode (candidate capacity at the 1800/s³ cap, `results/test_v3/REPORT.md`) is diagnosable and
   GT-free detectable. This is a real but secondary contribution, and reviewers will weigh it against a 0.10 F1 gap.
2. *Lower median volume error on EF grains* (no 1-voxel dilation). This is narrow and probably closed by a trivial
   calibration of Cellpose.
3. *Licensing / provenance.* Cellpose models are trained on CC-BY-NC data. For some industrial users that matters,
   but it is not a scientific argument.

**Defensible reframings (decision for the authors):**
- **(A) Benchmark / evaluation paper.** Make the contribution the *protocol and dataset*: 64 GT-free-defined crops
  over six specimens, the face-exclusion and ignore-rule analysis (the strict vs all-GT modes show how much the
  protocol choice moves the ranking), shape/contact fidelity, robustness. Report BA-TOS v3 as the
  classical/morphological reference and Cellpose as the learned reference.
- **(B) Hybrid.** Use the strengths of each: Cellpose for instances, connected operators for a measurable, auditable
  post-step (e.g. boundary calibration from the ToS, contact-area estimation, or GT-free QC that flags the crops where
  either method fails). Phase 3/5 (BA-TOS-prompted promptable models) fits here. It has to show a gain over cpsam_v2
  alone, which leaves at most ~0.002 F1 at IoU 0.50, so any gain would have to come in boundary or shape metrics.
- **(C) Harder regime.** The present data may simply be easy for a generalist model: well-separated, high-contrast
  grains, ~100 µm voxels, ~30 px diameter (Cellpose's training mean). §7 already shows that the obvious "harder"
  axes (noise, resolution) hurt **BA-TOS, not Cellpose**. v3 breaks at about 1.5–2× the native noise; cpsam_v2 does
  not move up to 10 %. A harder regime that favours BA-TOS would have to be found elsewhere: non-convex or fractured
  grains, materials or contrasts outside Cellpose's experience, or very large grains (> 120 px). Without such data,
  (C) is speculation, and the noise result argues against it.
- Do **not** reframe around the strict mode: it is explained above and would not survive review.

**Before any of this is written up:**
- Confirm the Cellpose-SAM dataset list in the preprint's methods (and ask the authors about the June 2026 models).
- Establish how the Pinzón GT was produced.
- ~~Finish cpsam / cpdino / cpdino-vitb.~~ Done: all four agree (§3b).
- If (B) is pursued, repeat §5 (shape/boundary) with cpdino, the Cellpose model with the best matched IoU, since any
  hybrid has to beat the strongest boundary baseline, not the default one.
- ~~Run the planned flow3D_smooth / diameter sensitivity (dev only).~~ Done: flat (§3c).
