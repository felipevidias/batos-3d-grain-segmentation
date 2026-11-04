# Recovered scripts inventory (linux_patches)

Written 2026-11-04 (Phase 0/1a, read-only).
Source: `C:\BA_TOS_DATA\linux_import\linux_patches\`. It was extracted from
`C:\Users\USER\Downloads\batos_linux_patches_v2.zip`, and all 94 files match the zip by SHA-256.
Identical copies of most bundles also exist under `C:\Users\USER\Desktop\batos codigos\`. That folder
also holds the Linux home dotfiles, `lysandre_code/experiments/crops/crop_600_800_300_200.tif` and
`lista_EFR_labels_py7zr.txt`.

None of these scripts can run as-is. They all assume the Linux layout
(`~/batos-3d-grain-segmentation/audit/`, `$REPO/.venv/bin/python`,
`/home/felipe/Pesquisa-Grain_Seg/Filtro_Conexo/ToSConOpCpp/CodeReport_LysandreMacke/...`). Every
pipeline runner also calls at least one **missing** script (see §0).

Conventions in this document:
- **In**: inputs (path patterns, dtype/layout). **Out**: outputs.
- Volumes are 3D TIFFs in (Z, Y, X) order, read with `np.squeeze(tifffile.imread(...))`. A trailing RGB axis is dropped where noted.
- "GT" = reference label volume (int labels, 0 = background), named `*_gt_clean.tif` in the E-series crop folders.

---

## 0. Missing scripts referenced by the recovered code

| Missing script | Where it lived (Linux) | Called by | Contract inferred from the callers |
|---|---|---|---|
| `batos_tos_seed_true_dual.py` | `audit/` | run_e11_all_crops, run_e12/e13_development_sensitivity, run_e13_e11_all_crops | **Core ToS candidate selector (Higra).** Positional: selector image (`*_flatten_lambda1p0_8u.tif`, uint8). Flags: `--out-dir --height-factor 0.012 --area-factor 22.0 --area-min 2 --area-max 2200 --delta-parent-min 0 --dark-percentile-max 78.0 --support-area-min 20 --support-area-max 45000 --support-max-steps 16 --support-area-ratio-preferred 18.0 --support-area-ratio-width 80.0 --support-contrast-norm 80.0 --delta-norm 20.0 --area-preferred-position 0.42 --area-preferred-width 0.55 --w-height 0.22 --w-delta 0.12 --w-area 0.18 --w-dark 0.18 --w-support 0.3 --min-score 0.35 --cleanup-open-radius 0 --cleanup-close-radius 0 --component-min 2 --component-max 1500 --max-selected-nodes 3600`. Writes `*_tosconop_autodual_candidates.csv` (one row per selected node; at least `node`, `polarity` ∈ {min,max}, `score`). The paper's "K = 1800 MIN + 1800 MAX" split must be inside this script. Not verifiable. |
| `e10_spatialize_candidates_no_gt.py` | `audit/` | all E11/E12/E13 runners | `--image <selector> --candidate-csv <above> --internal-source <internal_dark mask> --out-dir`. Writes `candidate_diagnostics_no_gt.csv` with columns `node, polarity, score, point_z, point_y, point_x, has_reconstructed_region, region_voxels, region_touches_volume_face, region_internal_overlap_fraction, region_dominant_internal_purity, region_internal_components`. Per paper §II-E, `region_voxels` is the **exclusive support** (each leaf voxel goes to its nearest selected ancestor). The gate reads it. |
| `audit_true_batos_e07.py` | `audit/` | run_e11_all_crops, run_e13_e11_all_crops | GT capacity audit (oracle coverage for MIN/MAX/DUAL at K). `--image --gt --candidate-csv --internal-source --out-dir --min-gt-area 1000`. Writes `variant_summary.csv` (including the `gt_with_{zero,one,two,three_or_more}_candidates` columns) and `variants/*_marker_points_uint16.tif`. Uses GT, so it is **diagnostic only** and not part of the pipeline. |
| `e09b_ignore_aware_audit.py` | `audit/` | every evaluation step | **Evaluator (seed-ignore protocol).** `--e09-dir <watershed dir with segmentations/*_distance_c1.tif> --gt --marker LABEL=normalized_marker.tif (repeat) --tos-label --out-dir`. Writes `evaluation/ignore_aware_results.csv` with `evaluation_protocol` ∈ {`seed_ignore`, `strict`, …}, `marker_source`, `{tp,fp,fn,precision,recall,f1}_iou{025,050,075}`, `valid_gt_count`. The seed-ignore rule per paper §III-C: a prediction whose **seed** lies inside an excluded GT is ignored, and a background-seeded prediction counts as an FP. Whether it removes predictions < 500 voxels is **unknown**. |
| `e10_*` (E10 frozen-holdout builder) | `audit/e10_frozen_holdout_*/` | run_e11_holdout_suite (inputs only) | It produced `frozen_protocol_and_inputs.json`, `e10_final_holdout_summary.json`, `02_spatial_no_gt/`, `03_frozen_markers/frozen_markers_uint16.tif` and `04_frozen_watershed_no_gt/segmentations/{tos_e10,baseline_existing}_distance_c1.tif`. Lost. |
| `batos_preprocess_opening_reconstruction_v2.py` | `CodeReport_LysandreMacke/code/segmentation/batos/` | e13_reconstruct_official_pipeline | `<gray crop> --out-dir --prefix --radius R` → `*_flatten_openrec_r{R}_8u.tif` (uint8). The algorithm is unknown. Best guess: background B = grayscale opening-by-reconstruction with a ball of radius R, then I − B + median(B), rescaled to uint8. |
| `run_batos_pipeline.sh` (Linux version) | `code/segmentation/batos/` | e13_reconstruct_official_pipeline | `<flat> <outdir> <lambda> <min_size>`. A same-named script exists in this repo (`scripts/run_batos_pipeline.sh`), but whether it is byte-identical is unverified. The repo copy calls 3 scripts that are also missing: `batos_tos_seed.py`, `filter_predicted_labels_by_size.py` and `count_tiff_labels.py` (that last one is in `src/batos/`). |
| `batos_tos_seed.py` | `code/segmentation/batos/` | run_batos_pipeline.sh step 2 | Older ToS seed builder (`--max-selected-nodes 900`, `--min-score 0.35`). Writes `*_tosconop_autodual_seed_raw.tif`. |
| `filter_predicted_labels_by_size.py` | same | run_batos_pipeline.sh step 6 | `--labels --out-dir --prefix --min-size --relabel` → `*_filtered_min{N}_labels.tif`. Trivial to rewrite. |
| `batos_make_paired_crop.py` | same | scripts/run_second_scan_…sh | `--gray --label --out-dir --x --y --z --size --prefix` → `*_gray_crop_8u.tif`, `*_label_crop.tif`. **The axis convention is unknown.** |
| `baseline_otsu_distance_watershed.py`, `run_macke_mintree_one_crop.sh` | same | scripts/run_second_scan_…sh | Legacy Otsu+WS and Macke Min-tree baselines (README, May 2026). |
| `batos_instance_metrics.py` | `code/segmentation/batos/` | legacy runner | **Present** as `src/batos/batos_instance_metrics.py` (Hungarian matching, `--ignore-gt-border`, `--min-gt-area`). It is not the E-series evaluator. |

Historical data directories that are also lost: `results_batos/` (12 crops, `yukiko_cleaned_batos_normal/crops/`, `e13_official_12crops/`), and every `audit/e1x_*` output root.

---

## 1. BA_TOS_E11_pre_article_tests (E11 + E12)

Purpose: regression guard against E10, source ablation, marker quality, reproducibility, 12-crop run, and E12 sensitivity and boundary audits.
The frozen settings stated in the README are K=3600, MIN, τs=0.35, NMS=6, EDT weight 0.20, MIN bonus 0.05 and connectivity 1.

| Script | What it does | CLI | In | Out | Depends on missing |
|---|---|---|---|---|---|
| `apply_e11_patches.sh` | Backs up, then copies `audit/*` into the repo | `<repo>` | — | files in `audit/` | — |
| `e11_build_markers_ablation.py` | **ToS marker builder (frozen).** Reads spatial diagnostics. Keeps rows with `has_reconstructed_region` and in-bounds points. Selects by polarity. Keeps candidates whose point falls in a 6-connected component of `internal_dark`. Priority = minmax(score) + 0.20·log1p(max EDT of that component)/log1p(global max EDT) + 0.05·[polarity==min]. Keeps **one candidate per internal component** (highest priority). Then greedy NMS with radius r (Euclidean, strict `<` r², ties by source row). Markers are **single voxels** at the candidate point. | `--diagnostic --internal-source --out-dir --polarity {min,max,dual} --nms-radius 6 --edt-priority-weight 0.20 --min-priority-bonus 0.05 [--expected-marker]` | `candidate_diagnostics_no_gt.csv`; `*_internal_dark_candidates.tif` (uint8 0/255) | `markers.csv`, `markers_uint16.tif` (point labels, zlib), `marker_config.json` | e10_spatialize (input CSV) |
| `e11_build_control_markers.py` | Non-ToS controls. `otsu_edt`: image → p1–p99 normalized → Otsu `>` thr → 6-conn closing ×1 → fill holes → drop components <500 → EDT → `peak_local_max(min_distance=1)` → NMS r → top-N by EDT. `internal_medial`: one EDT max per internal component. With `--match-count-from`, the count equals the number of 6-connected components in the reference marker. | `--image --mode {otsu_edt,internal_medial} --out-dir --nms-radius 6 [--internal-source] [--match-count-from \| --target-count]` | flat image (`*_flatten_openrec_r*_8u.tif`); internal mask | `markers_uint16.tif`, `markers.csv`, `control_config.json` | — |
| `e11_run_fixed_mask_watershed.py` | **Frozen watershed.** Mask = clean Otsu mask as above (no offset, 6-conn closing ×1, fill holes, components ≥500), OR the voxels of every `--mask-marker`. Relief = −(EDT normalized to p1–p99 inside the mask). Markers are clipped to the mask and relabeled with 6-connectivity. Then `skimage.watershed(connectivity=1, watershed_line=False)`. No size filter on predictions. | `--image --marker LABEL=TIFF (repeat) --mask-marker TIFF (repeat) --outside-policy {drop,error} --out-dir` | flat image (openrec); marker TIFFs | `frozen_mask_uint8.tif`, `negative_edt_float32.tif`, `normalized_markers/<label>.tif`, `segmentations/<label>_distance_c1.tif` (int32), `watershed_results.csv`, `fixed_mask_manifest.json` (SHA-256 of arrays) | — |
| `e11_marker_quality.py` | GT-based marker quality: coverage, exactly-one, marker precision, centrality | `--gt --marker LABEL=TIFF --out-dir --min-gt-area 1000` | GT + markers | `marker_quality_summary.{csv,json}`, `marker_per_gt.csv`, `marker_assignments.csv` | — |
| `e11_compare_regression.py` | Checks E10 vs E11 segmentation arrays and metrics (tol 1e-12) | `--e10-summary --e11-results --e10-tos-seg --e11-tos-seg --e10-baseline-seg --e11-baseline-seg [--tos-label --baseline-label --tolerance] --out` | E10 outputs | `regression_guard.json` | E10 artefacts |
| `e11_aggregate_all_crops.py` | Per-crop table, bootstrap CIs, paired deltas, Wilcoxon+Holm | `--root --out-dir --tos-label tos_e11_min --protocol seed_ignore --bootstrap 20000 --exclude-crop` | `*/evaluation/ignore_aware_results.csv`, `*/marker_quality/*` | `all_crop_results.csv`, `confirmatory_crop_results.csv`, `method_summary_bootstrap.csv`, `paired_source_ablation.csv`, `aggregate_summary.json` | e09b outputs |
| `e11_aggregate_capacity.py` | Aggregates the K=3600 oracle capacity | `--root --out-dir` | `*/capacity_k3600/variant_summary.csv` | `all_crop_capacity_k3600.csv`, `capacity_k3600_summary.csv`, `capacity_k3600_decision.json` | audit_true_batos_e07 |
| `e12_boundary_margin_audit.py` | **Strict protocol** evaluator. Valid GT = area ≥ min and bbox more than m voxels from every face. Ignore domain = invalid-GT voxels. Predictions are counted if they have any voxel in the domain. **Maximum-cardinality bipartite matching** on IoU ≥ thr edges. F1@0.50 and 0.75 by default. | `--gt --segmentation LABEL=TIFF --out-dir --min-gt-area 1000 --margins 0 2 4 8 12` | GT, segmentations | `boundary_margin_results.csv`, `boundary_margin_f1_pivot.csv`, `boundary_margin_summary.json` | — |
| `e12_physical_threshold_report.py` | Converts voxels to volume and d_eq | `--voxel-size-um --thresholds 500 1000 --out-dir` | — | `physical_thresholds.{csv,json}` | — |
| `run_e11_smoke_test.sh` | Synthetic dependency test | — | — | prints `E11 SMOKE TEST PASSED` | — |
| `run_e11_holdout_suite.sh` | E10 → E11 regression on the holdout and the full ablation | `[E10_DIR]` | `frozen_protocol_and_inputs.json`; `batos_openrec_lambda1_min500/04_markers/*` | `audit/e11_holdout_suite_*/` | E10 dir, e09b |
| `run_e11_reproducibility.sh` | Runs the holdout suite twice and compares | — | — | `reproducibility_report.json` | as above |
| `run_e11_all_crops.sh` | **Full 12-crop E11 pipeline.** Per crop it takes `01_preprocess_openrec_auto_minparam/*_flatten_openrec_r*_8u.tif` (FLAT, used for watershed and controls) and `batos_openrec_lambda1_min500/{01_preprocess/*_flatten_lambda1p0_8u.tif` (SELECTOR, used for the ToS), `03_background/*_internal_dark_candidates.tif`, `04_markers/*_batos_marker_components_raw.tif}`, plus `00_input/*_gt_clean.tif`. It then runs selector K3600 → capacity audit → spatialize → markers MIN/MAX/DUAL → edt_otsu and internal_medial controls → fixed-mask WS (mask markers = ToS MIN + baseline) → e09b → marker quality, and aggregates excluding crop04-dev. | `<CROPS_ROOT> [OUTROOT]` | see left | `audit/e11_all_crops_frozen_*/` | true_dual, e10_spatialize, e07 audit, e09b |
| `run_e12_development_sensitivity.sh` | τs/K/NMS sweep on `EFRGP01_00_crop04` | — (reads E10 protocol) | as above | `development_sensitivity_*.csv` | same |
| `run_e12_boundary_margin_all_crops.sh` | Runs E12 margins on every E11 crop | `<E11_ROOT> <CROPS_ROOT>` | `watershed/segmentations/*_distance_c1.tif` | `aggregate/boundary_margin_*.csv` | — |

## 2. BA_TOS_hotfix_v3

`apply_batos_hotfix_v3.sh <repo>` applies two shell fixes. (1) Loops used `PATH` as a loop variable, which broke `basename`. They now use `MARKER_PATH`. (2) The development crop is treated as an ID, not a path. No algorithmic change.

## 3. BA_TOS_E13_E16_conditional_acceptance (E13–E17)

| Script | What it does | CLI | In | Out | Missing deps |
|---|---|---|---|---|---|
| `apply_e13_e16.sh` | Installer, backup and Higra smoke test | `<repo>` | — | `audit/*` | — |
| `e13_build_12crop_manifest.py` | Symlinks the raw gray and GT of 12 historical crops into a normalized root. `CROP_MAP` gives the 12 IDs and origins (§2 of CLAUDE.md). | `--results-root --normalized-root [--replace-symlinks]` | `results_batos/EFRGP01_00_multicrop/crop0{2..5}_x…`, `EFRGP01_00_extra_crops/crop0{6,7}_…`, `EFRGP01_01_selected_crops/…`, each `00_input/*_gray_crop.tif\|*_gray_clean.tif` + `*_label_crop.tif\|*_gt_clean.tif` | `e13_12crop_manifest.csv` (scan, crop_id, gray, gt, crop_dir, source_crop_dir) | results_batos (data) |
| `e13_estimate_scale.py` | **GT-free scale estimate.** p1–p99 → uint8 → Otsu `>=` → closing ball(1) → fill holes → EDT → 5³ local maxima ≥ max(2, P75) → R = median of the top 1500 peaks (fallback P95, min 3). Derives `openrec_radius = max(8, round(1.5R))`, `lambda_after_openrec = 1` and `min_predicted_size = max(500, round(0.1·4/3πR³))`. | `--manifest --out-dir` | manifest | `<crop>_auto_params.json`, `e13_auto_scale_params_per_crop.csv`, `e13_auto_scale_global_median.json` | — |
| `e13_reconstruct_official_pipeline.py` | For each crop: (1) runs `batos_preprocess_opening_reconstruction_v2.py <gray> --radius R` → `*_flatten_openrec_r{R}_8u.tif`; (2) runs `run_batos_pipeline.sh <openrec> <tmp> 1 <min_size>`, so the legacy Gaussian flattening runs with **λ=1 on the openrec image**. Validates shapes, then renames atomically to `01_preprocess_openrec_auto_minparam/` and `batos_openrec_auto_minparam/`. | `--manifest --params-dir --project-root [--force] --expected-count 12` | manifest + params | per-crop folders; `e13_reconstructed_manifest.csv` with SHA-256 | openrec_v2, run_batos_pipeline.sh (+ batos_tos_seed.py, filter…) |
| `e13_validate_reconstruction.py` | **Fails if any selector is `lambda15p0`** or not `flatten_lambda1p0`, if the openrec name is unexpected, if a mask or markers are empty, or if shapes mismatch | `--manifest --out-dir` | reconstructed manifest | `e13_reconstruction_validation.{csv,json}`; prints `E13 RECONSTRUCTION VALIDATION PASSED: 12 CROPS` | — |
| `e13_aggregate_stratified.py` | Dev / holdout / retrospective strata and paired deltas | `--root --out-dir --bootstrap 20000 --tos-label` | `ignore_aware_results.csv` | `e13_stratified_*.csv/json` | e09b |
| `e14_build_mintree_candidates.py` | Min-tree candidate-source control capped at K=3600 (height·darkness·log-area score). It goes through the same spatialize/consolidation/NMS/WS steps. **This is not the Macke baseline.** | `--image --internal-source --out-dir --max-candidates 3600 --area-min 20 --area-max 45000` | selector image + internal | `candidate_diagnostics_no_gt.csv`, `mintree_source_control.json` | higra |
| `e15_make_crop_size_suite.py` | Central 128/160/200 sub-crops of the dev and holdout crops, keeping the parent scale params | `--e13-manifest --e13-params-dir --suite-root --sizes 128 160 200 [--force]` | E13 manifest + params | sub-crop gray/GT + `e15_crop_size_manifest.csv` | — |
| `e15_summarize_crop_size.py` | Crop-size tables | `--e11-root --out-dir` | `ignore_aware_results.csv` | `e15_crop_size_*.csv`, `e15_tos_crop_size_deltas.csv` | e09b |
| `e16_reviewer_gate_report.py` | Evidence matrix versus the reviews | `--e13-e11-root --sensitivity-root --boundary-summary --e15-summary-root --reproducibility-report [--physical-threshold-report] --out-dir` | aggregates | `e16_reviewer_gate_report.{md,csv,json}` | — |
| `run_e13_reconstruct_12crops.sh` | manifest → scale → reconstruct → validate | `[PROJECT_ROOT] [RESULTS_ROOT] [NORMALIZED_ROOT]` | | | as above |
| `run_e13_e11_all_crops.sh` | Same as `run_e11_all_crops.sh` but on the E13 root (`batos_openrec_auto_minparam/`) and with a `mintree_e14` source added | `<E13_ROOT> <E11_ROOT>` | | | true_dual, e10_spatialize, e07, e09b |
| `run_e13_development_sensitivity.sh` | **Source of Table IV.** On crop04-dev, sweeps τs ∈ {.25,.30,.35,.40,.45} (K3600, n6), K ∈ {900,1800,7200} (τ .35, n6) and NMS ∈ {4,8,10}. All variants share one fixed mask built from the default marker and the baseline. | `<E13_ROOT> [OUTROOT]` | | `development_sensitivity_{results,results_with_delta,compact}.csv` | true_dual, e10_spatialize, e09b |
| `run_e15_crop_size_sensitivity.sh` | Reconstruct + E11 on the size suite (**source of Table V**) | 4 roots | | | as above |
| `run_e13_full_pre_mac_suite.sh`, `run_e17_prepare_mac_handoff.sh` | Orchestration and handoff tarball | | | | |
| `README_EXECUTION.md`, `REVIEWER_RESPONSE_MATRIX.md` | Protocol and interpretation rules | | | | |

## 4. E23b–E23f (support-collapse gate, i.e. BA-TOS_sel)

| Script | What it does | CLI | In | Out | Missing deps |
|---|---|---|---|---|---|
| `e23b_degeneracy_gate_audit.py` | **Gate definition.** For each crop it merges `markers_tos_min/markers.csv` (node) with `spatial_no_gt/candidate_diagnostics_no_gt.csv` on `node` (one-to-one) and takes `region_voxels` of the **selected** nodes. g = median ≤ 3 ∧ q75 ≤ 4, using pandas `quantile` (linear interpolation). It then evaluates routes retrospectively from the existing `seed_ignore` F1 values. | `--e11-root --out-dir` | E11 root | `e23b_gate_by_crop.csv`, `e23b_triggered_crops.csv`, `e23b_route_aggregate.csv`, `e23b_summary.json`, `REPORT_E23B.md` | e10_spatialize outputs |
| `e23c_build_selective_marker_route.py` | **Router.** Same gate. If true it copies `markers_edt_otsu/markers_uint16.tif`, otherwise `markers_tos_min/markers_uint16.tif`, to `batos_selective_markers_uint16.tif` | `--candidate-diagnostics --tos-marker-csv --tos-marker-tiff --edt-marker-tiff --crop-id --out-dir` | as named | `route_manifest.json`, `route_summary.csv` | — |
| `run_e23c_selective_router_all12.sh` | Route → `e11_run_fixed_mask_watershed.py` (image = E13 openrec flat; mask markers = ToS MIN + **normalized** baseline) → e09b. Checks that the selected-source metrics and hashes equal the frozen E11 ones. README: mean F1@0.50 0.823817 → 0.906272, SD 0.367017 → 0.209900. | — (reads `audit/e13_stage2_last_root.txt`) | E11 and E13 roots | `audit/e23c_selective_router_all12_*/`, `e23c_selective_summary_all12.csv`, `e23c_decision.json` | e09b |
| `e23d_gate_robustness_audit.py` | Margin, bootstrap, threshold grid (2230 equivalent pairs) and scan-pair consistency. Writes the frozen config. | `--e11-root --e23c-summary --out-dir --bootstrap-replicates --seed` | | `E23D_FROZEN_GATE_CONFIG.json`, `e23d_*.csv`, `e23d_decision.json` | — |
| `e23e_finite_sample_stability.py` | LOO, subsampling at 90/80/70 %, adversarial flips (min fraction 0.113208) | `--e11-root --out-dir --replicates --seed` | | `e23e_*.csv`, `e23e_decision.json` | — |
| `e23f_build_frozen_handoff.py` | Hash-manifested code and report bundle | `--repo --out-dir` | | `HANDOFF_MANIFEST.json`, tarball | lists e09b, e10_spatialize as required |
| `install_*.sh`, `run_*.sh` | Copy into `audit/` and run with tee to a log | — | | | |

## 5. E24/E25 (reviewer closure)

| Script | What it does | CLI | Notes |
|---|---|---|---|
| `collect_e24_e25_pipeline_contract.py` | Read-only collector. Gathers script `--help` output and grep hits for τs/K/NMS/halo, and lists the E11 results. | `--repo --out-dir` | Its output tarball (`audit/e24_e25_pipeline_contract_*.tar.gz`) would contain the `--help` text of the missing scripts, but it is lost with the Linux machine. |
| `e24_reviewer_sensitivity_summary.py`, `run_e24_reviewer_parameter_sensitivity.sh` | Re-runs `run_e13_development_sensitivity.sh` and builds the reviewer table (**Table IV**) | `--sensitivity-root --out-dir` | |
| `e25_boundary_margin_readonly.py` | E12 margins on the frozen E11 segmentations | `--repo --e11-root --e13-root --out-dir` | Strict audit: 0.692 → 0.523 (paper) |
| `e25_reviewer_crop_size_summary.py` (3 identical copies), `run_e25_*` (original / fixed / resume / finalize) | Crop-size suite and tables (**Table V**) | `--boundary-root --e15-suite-root --e15-e11-root --out-dir` | The "fixed" run sets LC_ALL=C and adds disk checks. "finalize" infers `valid_gt_count` from the multiplicity columns. |
| `e11_aggregate_capacity_fixed.py` | Same as the E11 capacity aggregation plus that inference | `--root --out-dir` | |

## 6. batos_sibgrapi_v4_with_figures

An earlier draft (`main_batos_sibgrapi_2026_v4_with_figures.tex`) plus figures. It is superseded by the camera-ready source
`Downloads\2026_conf_sibgrapi_main_track_batos_felipe_dias_camera_ready\main.tex` (2026-09-15).

## 7. What exists in this repo (`src/batos`, `scripts`)

These are the legacy (May 2026, README-era) pipeline modules, and each is still needed by the E-series:
- `batos_preprocess.py`: Gaussian flattening per z-slice, σ = λ_xy, then I − B + median(B>0) and p0.5–p99.5 → uint8. Writes `*_flatten_lambda{λ}_8u.tif`. With λ=1 on the openrec image this produces the **E-series selector image**.
- `batos_background.py`: dark threshold = min(Otsu+offset, P_dark, max). Border-seeded 6-conn propagation gives `external_dark`, and `internal_dark = dark ∧ ¬external`. **This is the operational boundary rejection used by every E-series step** (ToS points must fall in `internal_dark` components). The pipeline uses `--dark-percentile 40 --otsu-offset -10 --dark-threshold-max 135`. The `--seed` input only affects its own candidate report, not `internal_dark_candidates.tif`.
- `batos_markers.py`: legacy internal-component markers (NMS 7, ≤420). These are the E-series `baseline_existing` markers, which also enlarge the frozen watershed mask.
- `batos_watershed_distance.py`: legacy watershed (Otsu − 10 on the robust uint8 image, closing ball(1), marker forcing). This is **not** the E-series watershed.
- `batos_instance_metrics.py`: legacy Hungarian matcher (ignores a prediction if its overlap with ignored GT exceeds its overlap with valid GT). It is **not** the seed-ignore protocol.
- `scripts/run_batos_pipeline.sh` has a cosmetic bug at line 33 (a literal `\n` inside `echo`) and calls 3 missing scripts.

---

## 8. Q5: `src/batos` on `main` vs `origin/batos-ls-auto-tile-pilot` (2026-11-04)

Method: `git show main:src/batos/<f>.py` vs `git show origin/batos-ls-auto-tile-pilot:CodeReport_LysandreMacke/code/segmentation/batos/<f>.py`, then `diff -w`. Nothing was checked out.

| Module | Result |
|---|---|
| `batos_preprocess.py`, `batos_watershed_distance.py`, `batos_instance_metrics.py`, `count_tiff_labels.py` | identical (ignoring whitespace) |
| `batos_background.py` | algorithmically identical. `main` lost the shebang and the 19-line docstring (`filter_tos_seed_background_reject_v14.py`, "Patch v14", which explains the border-connected background rejection). |
| `batos_markers.py` | algorithmically identical. `main` lost the shebang and the docstring (`internal_dark_to_markers_batos.py`; it notes the flattened input was "lambda=15"). |
| `batos_evaluate_against_gt.py` | **`main` is broken**: commit `324faaa` "Clean up imports" removed `argparse, csv, Path, numpy, tifffile`, so the script fails with NameError. Logic otherwise identical. Not fixed (frozen-code rule). New work must use the pilot-branch version or a new evaluator. |

Conclusion: the repo `src/batos` modules equal the Linux `CodeReport_LysandreMacke` versions on that branch. So `internal_dark` (boundary rejection) and `baseline_existing` markers can be regenerated from them. Caveat: the branch snapshot is the BA-TOS-LS-AUTO tile pilot commit `413920f`, and whether E13 used exactly this revision cannot be proven.

### Additional recovered files on `origin/batos-ls-auto-tile-pilot` (read-only via `git show`)

`CodeReport_LysandreMacke/code/segmentation/batos/`:
- `batos_tos_seed.py`: **legacy ToS selector**, Higra 3D ToS, the same flag set as the missing `batos_tos_seed_true_dual.py`. Writes `*_tosconop_autodual_{candidates.csv,seed_raw.tif,seed_centroids.tif,seed_report.txt}`. This is the basis for the v2 selector port.
- `batos_make_paired_crop.py`: crop extractor. `tifffile.memmap` then `arr[z0:z1, y0:y1, x0:x1]` from the CLI `--x --y --z --size`. It writes the **raw-dtype** gray crop, a robust-u8 copy and the label crop.
- `filter_predicted_labels_by_size.py` (+ `_labels_only`), `baseline_otsu_distance_watershed.py`, `batos_marker_diagnostics.py`, `batos_oracle_gtmask_watershed.py`, `diagnose_instance_errors.py`, `make_*` table/figure scripts, `paper_scripts/*` (BA-TOS-LS local split, full-image tile pilot, revision metrics).
- Also present: a C++ ToS library (`include/`, `src/*.cpp`, `CMakeLists.txt`). It is not used by the Python pipeline.

`CodeReport_LysandreMacke/results_batos/paper_tables/`: historical per-crop CSVs. Useful as **external checks**, for example:
- `batos_selected_crops.csv` / `batos_instance_selected_crops_summary.csv`: EFRGP01_00 crop02–07 legacy `gt_valid` (≥1000 voxels, not touching a face) = **265 / 216 / 272 / 129 / 271 / 263** (sum 1416), plus legacy BA-TOS F1@0.50 per crop.
- `batos_ls_auto_params/auto_params_per_crop.csv`: per-crop `n_labels` for all 12 crops (00: 394, 289, 398, 183, 367, 388; 01: 405, 289, 397, 187, 373, 390). **These count PREDICTED labels** of the legacy min500 BA-TOS output (`estimate_batos_ls_auto_params.py` reads `pred_path`), not GT. They are not usable as a GT check. For GT-derived counts see the next item.
- `batos_local_split_auto/EFRGP01_0X/cropNN/metrics_iou_0p50/*_instance_metrics.csv`, column `gt_valid` (legacy `batos_instance_metrics.py --min-gt-area 1000 --ignore-gt-border`), for **all 12 crops**:
  EFRGP01_00 crop02–07 = 265, 216, 272, 129, 271, 263 (Σ 1416); EFRGP01_01 crop02–07 = 262, 213, 264, 128, 265, 260 (Σ 1392). **Total 2808**, which equals the paper's valid-GT pool (Table III FN 2808 for "w/o boundary rejection"; strict audit 234.0 × 12). These are the reference counts for the Q7 axis check and the crop sanity check.
- `batos_local_split_large_labels/local_split_manifest.csv`: historical crop paths `results_batos/EFRGP01_0X_*/cropNN_xX_yY_zZ_s200/00_input/*_label_crop.tif`.
