# BA-TOS baseline reproduction (SIBGRAPI 2026 → `batos-paper-v2`)

**Status 2026-11-04: Phase 1 closed. Table II is reproduced (all statistics within ±0.0005) by running the frozen
E-series code verbatim. Tag `batos-paper-v2-frozen`.**

This is the "First Action" document required by `docs/CLAUDE_CODE_BATOS_JOURNAL_PLAN.md`. It is a short
entry point. Details are in:
- `docs/RECOVERED_SCRIPTS_INVENTORY.md`: every recovered script, its CLI, inputs and outputs.
- `docs/BATOS_PAPER_V2_SPEC.md`: stage-by-stage spec, the λ discrepancy analysis, best-guess register BG1–BG11, and round-4 replacement table (§3e).
- `docs/EXPERIMENT_LOG.md`: what was done, when, and with what result.

## 1. Which pipeline is the frozen SIBGRAPI baseline

The camera-ready paper (`paper/batos_sibgrapi_2026_camera_ready.pdf`) reports two pipelines. The frozen
baseline reproduced here is the **E-series pipeline** (E10 → E11 → E13 → E23). It produced Table II
(crop-macro), Table IV (sensitivity), Table V (crop extent) and the gate statistics:

1. GT-free scale estimate R̂ (`e13_estimate_scale`) → 8-bit p1–p99 → grayscale opening-by-reconstruction, ball radius max(8, round(1.5R̂)) → I − B + median(B) → p1–p99 8-bit → Gaussian flattening λ = 1 (selector image).
2. Higra 3D Tree of Shapes, true-dual selector (`batos_tos_seed_true_dual.py`): 1800 MIN + 1800 MAX candidates, |height|, score ≥ τs = 0.35, hierarchy suppression.
3. `e10_spatialize`: exclusive supports via the nearest selected ancestor; boundary rejection via the internal dark component.
4. `e11_build_markers_ablation`: MIN markers, priority = score + 0.20·EDT + 0.05·MIN, NMS radius 6. EDT-Otsu control at matched count.
5. Support-collapse gate (`e23c`): median(v_i) ≤ 3 ∧ Q75(v_i) ≤ 4 → EDT-Otsu markers.
6. `e11_run_fixed_mask_watershed`: Otsu mask (no offset, closing ×1, fill holes, components ≥ 500, ∪ marker voxels) → watershed on −EDT, connectivity 1.
7. **No prediction size filter** in the published numbers (established in Phase 1 closure, §6).
8. `e09b` seed-ignore evaluation (paper §III-C). The frozen script evaluates IoU 0.50/0.75; 0.25 is added by `scripts/e09b_with_iou025.py`, which overrides only the default threshold tuple.

The May-2026 pipeline (λ_xy = 15, Otsu −10, pooled F1 0.974 in Table III) is a documented legacy preset
only (`batos-legacy-may2026`). Table I of the paper mixes the two (see CLAUDE.md §2).

## 2. Inputs (verified)

| Item | Location |
|---|---|
| Greyscale volumes | `C:\Users\USER\Downloads\PhDUpload\PhDUpload\EFR\Greyscale\EFRGP01_0{0,1}.tif` (uint16, Z,Y,X = 1600,1400,1400) |
| GT label volumes | `C:\BA_TOS_DATA\gt_full\EFRGP01_0{0,1}.tif` (from `EFR_labels.7z`; uint16). `Desktop\BATOS_DATA\*` is the same GT. |
| Crop definitions | `data/metadata/crops_sibgrapi.json` (12 crops, 200³, axis order `zyx` **verified**) |
| Cut crops | `C:\BA_TOS_DATA\crops\<crop_id>\{gray.tif, labels.tif}`; manifest `data/metadata/crops_manifest.csv` |
| Valid GT | ≥ 1000 voxels and no face contact; 2808 in total, equal to the historical per-crop counts |
| Frozen code | `C:\BA_TOS_DATA\linux_import\chat_recovered\e13_crop05_metadata_full_20260725_182408\audit\` and `zsh_recovered\`; pilot `batos_tos_seed.py` staged at `C:\BA_TOS_DATA\frozen_code\pilot_branch_413920f\`. SHA-256 of every script is recorded per step in `frozen_run.json`. |
| Environment | `C:\BA_TOS_DATA\envs\batos-py312` (Python 3.12.10; numpy 2.5.1, scipy 1.18.0, scikit-image 0.26.0, higra 0.6.13, pandas 3.0.5, tifffile 2026.7.14, i.e. the Linux pins), frozen in `environments/batos-py312.txt` |

## 3. Recovered vs reconstructed

Everything in the reproduction of record is **recovered frozen code run unmodified** (`src/batos_v2/frozen.py`
only orchestrates subprocesses with the flags and stage order of `e13_stage2_worker_*.sh` / `run_e13_e11_all_crops.sh`).
The rounds 1–3 port (`batos_v2.run_crop`, `preprocess.py`, `tos_select.py`, …) and its best guesses BG1–BG11
are kept for reference only; spec §3e lists which guess each recovered script replaced.

Not run: `audit_true_batos_e07.py` (GT capacity audit, not recovered), the internal_medial and E14 Min-tree
controls. Hence the Min-tree-core+WS baseline (paper 0.768) is **not compared** in Phase 1.

## 4. Fidelity targets and outcome

| Check | Target | Reproduced |
|---|---|---|
| EFRGP01_00_crop04 F1@0.50 / F1@0.75 | 0.985 / 0.978 | 0.9852 / 0.9778 |
| Gate | fires only on both crop05 with (2, 3); closest non-triggered (17, 43) | yes: both crop05 (2, 3); closest 01_crop03 (17, 43) |
| Crop-macro F1@0.50 ToS | 0.824 ± 0.367, median 0.985, min 0.012 | 0.8238 ± 0.3670, 0.9850, 0.0121 |
| Crop-macro F1@0.50 sel | 0.906 ± 0.210, median 0.985, min 0.260 | 0.9063 ± 0.2099, 0.9850, 0.2599 |
| Crop-macro F1@0.50 Otsu+WS | 0.830 ± 0.182, median 0.882, min 0.260 | 0.8301 ± 0.1819, 0.8816, 0.2599 |
| Crop-macro F1@0.75 ToS → sel | 0.817 → 0.899 | 0.8175 → 0.8992 |
| Retained ToS markers per crop | 53–405 | 53 (00_crop05) – 405 (00_crop02), all 12 crops |
| K / τs / NMS sweeps (Table IV), crop extent (Table V) | see CLAUDE.md §2 | **not run in Phase 1** |
| Pooled F1@0.50 | recomputed, reported next to legacy 0.974 | sel 0.9496, ToS 0.9116 (no filter); not a pass/fail target |

SDs are **sample** SDs (ddof = 1); this is the paper's convention (population SDs are 0.008–0.016 lower).
Results: `results/batos_sibgrapi/` (git-ignored; regenerate with the commands below).

## 5. Commands

Reproduce one crop (reproduction of record; refuses to overwrite an existing `frozen_run.json`):
```
C:\BA_TOS_DATA\envs\batos-py312\Scripts\python.exe -I -m batos_v2.frozen --crop-id EFRGP01_00_crop04 [--runs-root C:\BA_TOS_DATA\runs\batos_paper_v2_frozen]
```
Aggregate the 12 crops (Table II, pooled, size filter on/off, `results/all_runs.csv`):
```
C:\BA_TOS_DATA\envs\batos-py312\Scripts\python.exe -I scripts\aggregate_sibgrapi.py
```
Run both from the repo root. `batos_v2` is installed editable in the venv, so `-I` works. Tests:
`python -I -m unittest discover -s tests -t .`. Runtime is 5–8 min per crop (290–490 s) on a Ryzen 7 5700X (CPU only), and peak memory is 3.3–4.8 GB per crop,
set by the openrec stage. Runs are deterministic: 12 full reruns were byte-identical.

## 6. Prediction size filter (500 voxels)

Applied after the watershed with the recovered `filter_predicted_labels_by_size.py --min-size 500` **without
`--relabel`**. e09b links a prediction to its seed through the label value (watershed label = marker label), so sequential
relabelling corrupts the seed-ignore sets. A first aggregation with `--relabel` gave ToS crop04 0.726 with only 4 objects removed and was discarded.

| Route | Crop-macro F1@0.50, filter off | filter 500 | Pooled F1@0.50, filter off | filter 500 |
|---|---|---|---|---|
| BA-TOS_ToS | 0.8238 | 0.8297 | 0.9116 (P 0.932, R 0.892) | 0.9351 (P 0.982, R 0.892) |
| BA-TOS_sel | 0.9063 | 0.9088 | 0.9496 (P 0.963, R 0.937) | 0.9522 (P 0.968, R 0.937) |
| Otsu+WS | 0.8301 | 0.8301 | 0.8657 | 0.8657 (no object < 500) |

The filter removes only false positives (TP unchanged on every crop). Table II is matched by the **filter-off** run
(ToS min 0.012 vs 0.0153 with the filter), so the published E-series numbers did not use it.
