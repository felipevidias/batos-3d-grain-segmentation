# CLAUDE.md — BA-TOS Journal Extension

> This file lives at the root of `~/projects/batos-3d-grain-segmentation`.
> Claude Code reads it automatically at every session. The full experimental plan lives in
> `docs/CLAUDE_CODE_BATOS_JOURNAL_PLAN.md` — read that file before any non-trivial task.

## 1. What this project is

Repository: https://github.com/felipevidias/batos-3d-grain-segmentation
Local clone: `~/projects/batos-3d-grain-segmentation`
Local data (not in the repo): `~/Documents/BA_TOS_DATA/EFRGP01_00.tif`, `~/Documents/BA_TOS_DATA/EFRGP01_01.tif`

We are extending the SIBGRAPI 2026 paper
**"BA-TOS: Boundary-Aware Tree-of-Shapes Markers for 3D Grain Instance Segmentation in Micro-CT"**
(Dias, Patrocínio Jr., Passat, Kenmochi, Guimarães) into a journal submission.

BA-TOS = conservative 3D flattening → Tree-of-Shapes (ToS) dark MIN-like candidate extraction →
rejection of crop-boundary-connected candidates → internal marker refinement →
GT-free support-collapse safeguard (gate `median(v_i) ≤ 3 ∧ Q75(v_i) ≤ 4` → EDT-Otsu markers) →
marker-controlled distance watershed inside the Otsu mask → size filter (< 500 voxels removed).

Current repo layout (as published):

```
src/batos/        core code (e.g. batos_instance_metrics.py, count_tiff_labels.py)
scripts/          experiment scripts
sample_data/      two SYNTHETIC crops for smoke tests only — NOT the paper data
paper/figs/       figures
requirements.txt  numpy scipy scikit-image tifffile pillow matplotlib pandas py7zr higra
```

## 2. Frozen SIBGRAPI facts (never silently change these)

| Item | Value |
|---|---|
| Data | EFRGP01_00.tif, EFRGP01_01.tif (Zenodo 8014905), ~100 µm/voxel |
| Crops | 6 per scan × 2 scans = 12 crops of 200³. IDs are `EFRGP01_{00,01}_crop02..07` (the paper and the scripts both use this numbering). Origins (x,y,z), size 200: crop02 600/800/300, crop03 900/600/300, crop04 500/500/700, crop05 800/900/700, crop06 350/650/500, crop07 700/400/500 — identical for both scans. Development crop = EFRGP01_00_crop04; frozen temporal holdout = EFRGP01_01_crop04. |
| Flattening (`batos-paper-v2`, E-series) | opening-by-reconstruction, radius max(8, round(1.5·R̂)) with R̂ from the GT-free scale estimate (`e13_estimate_scale.py`), then Gaussian flattening λ = 1 → selector image. **Not** λxy = 15 (see the note below the table). |
| Flattening (legacy preset `batos-legacy-may2026` only) | λxy = 15 (as printed in Table I of the paper) |
| ToS area factor / range | 22 / 2–2200 |
| ToS score threshold τs | 0.35 |
| Candidate capacity K | 3600 (1800 MIN + 1800 MAX before eligibility) |
| NMS radius | 6 voxels |
| Boundary rule | reject C if C ∩ ∂Ω ≠ ∅ (any of the 6 faces) |
| Gate | m_v ≤ 3 and q_v ≤ 4 → EDT-Otsu fallback (fires only on `*_crop05` of both scans) |
| Mask (`batos-paper-v2`) | Otsu **without offset** on the p1–p99-normalized openrec image, 6-conn closing ×1, fill holes, components ≥ 500, ∪ ToS-MIN and baseline marker voxels (`e11_run_fixed_mask_watershed.py`) |
| Mask (legacy preset only) | Otsu, offset −10, closing 1 voxel (as printed in Table I) |
| Filters | pred < 500 voxels removed (**v2: config switch, unknown whether the published E-series numbers used it; Phase 1c runs both flat-500 and off**); GT < 1000 voxels or touching a face excluded |
| Evaluation (v2) | seed-ignore, paper §III-C: a prediction seeded exclusively inside an excluded reference component is ignored; a background-seeded prediction remains a false positive |
| IoU thresholds | 0.25 / 0.50 / 0.75 |
| Expected results (v2 fidelity targets) | crop-macro mean F1@0.50: BA-TOS_ToS 0.824 ± 0.367 (median 0.985, min 0.012); BA-TOS_sel 0.906 ± 0.210 (min 0.260); mean F1@0.75 0.817 → 0.899 |
| Legacy result (not a v2 target) | pooled core F1@0.50 = 0.974 (P 0.980, R 0.967, FP 39, FN 77; Table III and README of 2026-05-08). It came from the May-2026 λ = 15 pipeline. **The pooled score will be recomputed with v2 and reported next to 0.974, not forced to match it.** |
| Baselines in paper | Min-tree-core+WS (0.768 macro) = E14 Min-tree source control (`e14_build_mintree_candidates.py`); Otsu+WS (0.830 macro) = `edt_otsu` matched-count control (`e11_build_control_markers.py`). Decided 2026-11-04. |
| Sensitivity (crop04 of scan 00) | τs 0.25–0.45 flat (0.985); K = 900/1800/3600/7200 → 0.510/0.833/0.985/0.998; NMS 4/6/8/10 → 0.976/0.985/0.985/0.985 |
| Crop extent | 128³ / 160³ / 200³ → mean F1@0.50 0.851 / 0.967 / 0.985 |

Reproduction is only "done" when these numbers are matched (within the fidelity criteria in Section 4, Phase 1) and saved under `results/batos_sibgrapi/`.

> **Note on Table I of the paper (decided 2026-11-04).** Table I mixes two pipelines.
> "Flattening λxy = 15" and "Mask Otsu offset −10 / closing 1 voxel" describe the **legacy May-2026
> pipeline**, which produced only the pooled ablation (Table III, pooled 0.974). The crop-macro
> results (Table II), the sensitivity (Table IV), the crop extent (Table V) and the gate statistics come from the
> **E-series pipeline** (opening-by-reconstruction + λ 1, Otsu mask without offset, K = 3600 MIN,
> τs 0.35, NMS 6, exclusive-support gate, seed-ignore evaluation). That is what `batos-paper-v2`
> implements. Evidence: `docs/BATOS_PAPER_V2_SPEC.md` §1. The legacy pipeline is kept only as the
> documented preset `batos-legacy-may2026`.

### 2b. Where the camera-ready code actually lives (found 2026-11-05)

The SIBGRAPI revision was run on the Linux machine, as a series of patch packages applied to
`~/batos-3d-grain-segmentation/audit/` there (E10 frozen markers → E11 regression/ablation →
E12 sensitivity & boundary audit → E13 reconstruction of preprocessing on 12 crops → E14 Min-tree
control → E15 crop-size suite → E16 reviewer gate → E17 Mac handoff). Copies of E11, E12, E13–E16
and hotfix v3 are in `linux_patches/`. Confirmed frozen settings there: K=3600, MIN candidates,
τs=0.35, NMS=6, EDT weight 0.20, connectivity 1, seed-ignore evaluation at IoU 0.25/0.50/0.75.
Status 2026-11-05 (update): recovered. GT = `C:\Users\USER\Downloads\PhDUpload\PhDUpload\EFR\Label\EFR_labels.7z` (2.44 GB, 55 full-volume label TIFFs EFRGP01_00..54; crop GT = 200³ sub-blocks of EFRGP01_00/01 at the origins in §2). Gate = E23b/E23c (`e23b_degeneracy_gate_audit.py`, `e23c_build_selective_marker_route.py`: median(region_voxels)≤3 ∧ q75≤4 → EDT-Otsu); E23d robustness, E23e finite-sample, E23f frozen handoff, E24 sensitivity, E25 crop size — all in `linux_patches/`. The Linux research folders were deleted (`rm -rf ~/batos-3d-grain-segmentation`, `rm -rf ~/Pesquisa-Grain_Seg`), so still missing: the support-collapse gate (mv≤3 ∧ qv≤4 → EDT-Otsu) and the BA-TOS_ToS vs
BA-TOS_sel split (must be a package after E17); `audit/batos_tos_seed_true_dual.py` (core ToS
seed builder); `CodeReport_LysandreMacke/code/segmentation/batos/*` (preprocessing, pipeline,
metrics); `results_batos/` (GT labels, 12 crops, E13 outputs). Open discrepancy: E13 preprocessing
is opening-by-reconstruction + λ=1.0 and explicitly rejects `lambda15p0`, whereas the paper's
Table I says λxy=15 — resolve before Phase 1. **Resolved 2026-11-04:** v2 follows E13 (see the §2 note).
**Partial recovery 2026-11-04:** branch `origin/batos-ls-auto-tile-pilot` holds
`CodeReport_LysandreMacke/code/segmentation/batos/*` (incl. `batos_tos_seed.py` — the direct ancestor of
`batos_tos_seed_true_dual.py` with the identical flag set — `batos_make_paired_crop.py`,
`filter_predicted_labels_by_size.py`, `baseline_otsu_distance_watershed.py`) and historical per-crop tables
under `CodeReport_LysandreMacke/results_batos/paper_tables/`. Read with `git show`; never check out over `main`.
Still missing: `batos_tos_seed_true_dual.py`, `e10_spatialize_candidates_no_gt.py`, `e09b_ignore_aware_audit.py`,
`batos_preprocess_opening_reconstruction_v2.py`, `audit_true_batos_e07.py`.
Canonical data: `C:\Users\USER\Downloads\PhDUpload\PhDUpload\EFR\{Greyscale,Label}` (the `~/Documents` paths in §1 are
from the Mac and do not exist here). `Desktop\BATOS_DATA\*.tif` are **not greyscale**: they are byte-identical to the GT label volumes EFRGP01_00/01 from
`EFR_labels.7z` (extracted to `C:\BA_TOS_DATA\gt_full\`). Crops: axis order verified `zyx`
(`vol[z:z+200, y:y+200, x:x+200]`); cut to `C:\BA_TOS_DATA\crops\<crop_id>\{gray,labels}.tif`, manifest
`data/metadata/crops_manifest.csv` (valid GT 2808, matching the historical per-crop counts exactly).

## 3. Hard rules

1. **Never modify the frozen BA-TOS code or overwrite its outputs.** New work goes in new folders
   (`experiments/batos_v2`, `experiments/batos_fullvolume`, `experiments/batos_hybrid`, `baselines/<name>`).
2. Inspect before modifying; say which file you are changing and why.
3. Never invent paths, metrics, or results. If GT or metadata are missing, **fail loudly**.
4. Do not tune any parameter per crop using ground truth. Oracle-prompt runs are upper bounds and must be labelled `oracle prompt upper bound`.
5. Keep pooled-object and crop-macro metrics in separate, clearly labelled tables.
6. Every run writes: command, timestamp, git commit, environment, parameters, inputs, outputs, metrics → one machine-readable `results/all_runs.csv` (or parquet) with the columns in the plan.
7. Prefer YAML configs in `configs/` over hard-coded settings.
8. Keep a running `docs/EXPERIMENT_LOG.md`. Append, never rewrite.
9. Data lives locally in `~/Documents/BA_TOS_DATA/`; never run heavy jobs from the SMB share (`/Volumes/BATOS_DATA`).
10. Never run Cellpose or any foundation model on the full ~6 GB volume as a first experiment; start from the 12 crops.
11. `sample_data/` is synthetic and only for smoke-testing scripts.
12. Do not commit `.tif`, `.tiff`, `.7z`, `.zip`, model weights or result folders.

## 4. Roadmap (do in order; do not skip ahead)

| Phase | Goal | Done when |
|---|---|---|
| 0 | Inspect repo, write `docs/BA_TOS_BASELINE_REPRODUCTION.md` | Document lists entry points, crop definitions, GT paths, frozen params, deps, exact command for one crop |
| 1 | **Re-implement** BA-TOS as Algorithm 1, following the **E-series** pipeline rather than Table I's legacy entries (see the note in §2; the camera-ready code was lost with the Linux machine; see 2b). New preset `batos-paper-v2` in batos-framework; the E11 ablation scripts in `linux_patches/` are the interface reference. Openrec radius, λ and the prediction size filter are YAML-configurable; every reconstructed step is marked "best guess" in `docs/BATOS_PAPER_V2_SPEC.md`. Sub-phases: 1a spec, 1b environment + data (venv, GT extraction, 12 crops), 1c implementation + validation. | Fidelity criteria: EFRGP01_00_crop04 F1@0.50 ≈ 0.985; both crop05 trigger the gate with mv=2, qv=3 and the other ten do not (closest non-triggered: mv=17, qv=43); K sweep on crop04 ≈ 0.510/0.833/0.985/0.998; crop-macro mean F1@0.50 ≈ 0.906 (sel) / 0.824 (ToS). Reported for size filter = 500 **and** off. The pooled F1@0.50 is recomputed with v2 and reported next to the legacy 0.974; it is not a pass/fail target. Any deviation is reported, not hidden. Result frozen under `results/batos_sibgrapi/` and tagged `batos-paper-v2-frozen`. |
| 2 | Cellpose baselines (`cpsam_v2`, `cpdino`, `cpdino-vitb`) | Default settings, same 12 crops, same evaluator; then small fixed sensitivity (`flow3D_smooth` 0/1/2) |
| 3 | 3D promptable models: nnInteractive → SAM-Med3D(-turbo) → Medical-SAM3 | Three regimes recorded separately: automatic / oracle-prompt / BA-TOS-prompted |
| 4 | LocateAnything (2D grounding only) on slices z = 50, 75, 100, 125, 150; fixed prompts | Box P/R/IoU, missed grains, duplicates |
| 5 | Improve BA-TOS: adaptive K, continuous support-quality gate, tiled full-volume BA-TOS, hybrid BA-TOS→prompts→FM | Each as its own experiment folder; gates frozen on dev data before test |
| 6 | Journal-grade evaluation: P/R/F1@{.25,.5,.75}, AP50, AP75, mean matched IoU, over/under-seg, bootstrap CIs; morphological errors (volume, d_eq, centroid, surface, sphericity, elongation, orientation); runtime & memory table | `scripts/aggregate_results.py` and `scripts/morphology_analysis.py` produce all tables/figures from `results/all_runs.csv` |

Data splits for the journal: **Development** = the 12 SIBGRAPI crops; **Test** = new crops never used for design;
**External** = another scan / material / acquisition condition. Freeze parameters before touching Test.

## 5. Scripts to create (Phase 1–2)

```
scripts/run_batos_baseline.py   # frozen SIBGRAPI config, one or all 12 crops
scripts/run_cellpose.py         # one YAML config → label volume per crop
scripts/run_evaluation.py       # wraps src/batos/batos_instance_metrics.py, ignore-aware protocol
scripts/run_all_baselines.py
scripts/aggregate_results.py
scripts/morphology_analysis.py
configs/batos_sibgrapi.yaml
configs/cellpose_cpsam_v2.yaml
configs/cellpose_cpdino.yaml
```

## 6. References already cited in the SIBGRAPI paper

Keep this list; the journal version adds Cellpose-SAM, nnInteractive, SAM-Med3D, Medical-SAM3, LocateAnything, plus grain-morphology / micro-CT segmentation literature. Do not alter these entries.

1. G. Pinzón, E. Andò, A. Tengattini, G. Viggiani, "X-ray tomography analysis of particle morphology influence on granular deformation processes," *Open Geomechanics*, vol. 6, p. A2, 2025.
2. G. Pinzón, E. Andò, A. Tengattini, G. Viggiani, "4D image data and experimental results for granular deformation processes," Zenodo record 8014905, 2023. https://zenodo.org/records/8014905
3. L. Macke, "Analyse topologique de mouvements de grains dans une séquence d'images 3D," Master 2 internship report, Université Gustave Eiffel / GREYC, Université de Caen Normandie, 2024.
4. P. Salembier, A. Oliveras, L. Garrido, "Anti-extensive connected operators for image and sequence processing," *IEEE TIP*, vol. 7, no. 4, pp. 555–570, 1998.
5. E. Carlinet, T. Géraud, "A comparative review of component tree computation algorithms," *IEEE TIP*, vol. 23, no. 9, pp. 3885–3895, 2014.
6. S. Beucher, F. Meyer, "The morphological approach to segmentation: The watershed transformation," in *Mathematical Morphology in Image Processing*, E. R. Dougherty, Ed., Marcel Dekker, 1992, pp. 433–481.
7. C. Ballester, V. Caselles, P. Monasse, "The tree of shapes of an image," *ESAIM: COCV*, vol. 9, pp. 1–18, 2003.
8. N. Passat, Y. Kenmochi, "A topological tree of shapes," in *Discrete Geometry and Mathematical Morphology*, Springer, 2022, pp. 221–235.
9. J. Mendes Forte, N. Passat, A. Shimizu, Y. Kenmochi, "Consistent connected operators based on trees of shapes," *SIAM J. Imaging Sciences*, vol. 18, no. 4, pp. 2547–2579, 2025.
10. C. Stringer, T. Wang, M. Michaelos, M. Pachitariu, "Cellpose: A generalist algorithm for cellular segmentation," *Nature Methods*, vol. 18, pp. 100–106, 2021.
11. N. Otsu, "A threshold selection method from gray-level histograms," *IEEE Trans. SMC*, vol. 9, no. 1, pp. 62–66, 1979.
12. B. Perret, G. Chierchia, J. Cousty, S. J. F. Guimarães, Y. Kenmochi, L. Najman, "Higra: Hierarchical graph analysis," *SoftwareX*, vol. 10, p. 100335, 2019.

## 7. Kickoff prompt (paste as the first message in Claude Code)

```
Read CLAUDE.md and docs/CLAUDE_CODE_BATOS_JOURNAL_PLAN.md fully before doing anything.

Then do ONLY Phase 0:
1. Inspect this repository (src/batos, scripts, requirements.txt, README.md).
2. Check that ~/Documents/BA_TOS_DATA/EFRGP01_00.tif and EFRGP01_01.tif exist locally;
   report shape, dtype and size with tifffile WITHOUT loading the whole volume into RAM.
3. Report Mac hardware (chip, RAM, GPU/MPS availability) and the Python/Conda setup.
4. Identify the exact scripts for crop extraction, flattening, ToS construction (Higra),
   marker selection, gate, watershed, evaluation and metric computation.
5. Identify where the 12 crop definitions (origins) and the GT label files live, and their format.
6. Identify the frozen SIBGRAPI parameters in the code and confirm they match Section 2 of CLAUDE.md.
   If anything differs, STOP and report it; do not "fix" it.
7. Write docs/BA_TOS_BASELINE_REPRODUCTION.md with all of the above and the exact command
   to reproduce ONE crop (EFRGP01_00_crop04).
8. Create docs/EXPERIMENT_LOG.md and log this session.

Do not modify any existing code or results. Do not install Cellpose yet.
At the end, list the open questions you could not resolve from the repo alone.
```
