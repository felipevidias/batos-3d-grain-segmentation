# BA-TOS baseline reproduction (SIBGRAPI 2026 → `batos-paper-v2`)

This is the "First Action" document required by `docs/CLAUDE_CODE_BATOS_JOURNAL_PLAN.md`. It is a short
entry point. Details are in:
- `docs/RECOVERED_SCRIPTS_INVENTORY.md`: every recovered script, its CLI, inputs and outputs, and what is missing.
- `docs/BATOS_PAPER_V2_SPEC.md`: stage-by-stage spec, the λ discrepancy analysis, and the best-guess register BG1–BG11.
- `docs/EXPERIMENT_LOG.md`: what was done, when, and with what result.

## 1. Which pipeline is the frozen SIBGRAPI baseline

The camera-ready paper (`paper/batos_sibgrapi_2026_camera_ready.pdf`) reports two pipelines. The frozen
baseline reproduced here is the **E-series pipeline** (E10 → E11 → E13 → E23). It produced Table II
(crop-macro), Table IV (sensitivity), Table V (crop extent) and the gate statistics:

1. GT-free scale estimate R̂ → grayscale opening-by-reconstruction flattening, radius max(8, round(1.5R̂)) → Gaussian flattening λ = 1 (selector image).
2. Higra 3D Tree of Shapes, MIN-/MAX-like candidates, score ≥ τs = 0.35, capacity K = 3600 (1800 + 1800).
3. Boundary rejection: candidates must lie in an internal dark component (a dark region not connected to any of the 6 faces).
4. One point marker per internal component, priority = score + 0.20·EDT + 0.05·MIN, NMS radius 6.
5. Support-collapse gate: median(v_i) ≤ 3 ∧ Q75(v_i) ≤ 4 over exclusive supports → EDT-Otsu markers (matched count).
6. Fixed Otsu mask (no offset, closing ×1, fill holes, components ≥ 500, ∪ marker voxels) → watershed on −EDT, connectivity 1.
7. Prediction size filter 500 (whether it was applied is unknown, so it is run both on and off).
8. Seed-ignore evaluation at IoU 0.25 / 0.50 / 0.75 (paper §III-C).

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
| Environment | `C:\BA_TOS_DATA\envs\batos-py312` (Python 3.12.10, higra 0.6.13), frozen in `environments/batos-py312.txt` |

## 3. What was recovered vs reconstructed

**Recovered and ported as-is:** the E13 scale estimator, the E11 marker builder (consolidation + NMS), the E11 control
markers (EDT-Otsu), the E11 fixed-mask watershed, the E23b/c gate and router, and the E12 strict evaluator (matching). Also
`src/batos/batos_preprocess.py` (λ flattening), `batos_background.py` (internal dark mask), `batos_markers.py` (baseline
markers, used only to widen the mask), and `batos_tos_seed.py` on `origin/batos-ls-auto-tile-pilot` (ToS selector).

**Reconstructed (best guesses, all configurable in `configs/batos_paper_v2.yaml`):** BG1–BG11 in the spec §3b.
The main ones:
- BG3–BG5: the opening-by-reconstruction operator and how its output is used.
- BG8: the true-dual K split.
- BG9: exclusive supports.
- BG10: the size filter.
- BG11: the seed-ignore details.

The original `e10_spatialize`, `e09b` and `true_dual` scripts are lost.

## 4. Fidelity targets

| Check | Target |
|---|---|
| EFRGP01_00_crop04 (development) F1@0.50 / F1@0.75 | 0.985 / 0.978 |
| Gate | fires only on both crop05 with (m_v, q_v) = (2, 3); closest non-triggered (17, 43) |
| K sweep on crop04 (900/1800/3600/7200) | F1@0.50 0.510 / 0.833 / 0.985 / 0.998 |
| τs sweep 0.25–0.45 / NMS 4,6,8,10 on crop04 | flat 0.985 / 0.976, 0.985, 0.985, 0.985 |
| Crop-macro F1@0.50 (12 crops) | BA-TOS_ToS 0.824 ± 0.367; BA-TOS_sel 0.906 ± 0.210 |
| Crop-macro F1@0.75 | 0.817 → 0.899 |
| Retained ToS markers per crop | 53–405 |
| Pooled F1@0.50 | recomputed and reported next to legacy 0.974; not a pass/fail target |

Each target is reported with the size filter on (500) and off. Deviations are reported, never tuned away.

## 5. Command to reproduce one crop

```
C:\BA_TOS_DATA\envs\batos-py312\Scripts\python.exe -I -m batos_v2.run_crop --crop-id EFRGP01_00_crop04 --config configs\batos_paper_v2.yaml
```
Run it from the repo root. `batos_v2` is installed in editable mode in the venv (`pip install -e .`, see
`pyproject.toml`), so `-I` works. Tests: `python -I -m unittest discover -s tests -t .`. Outputs go to
`C:\BA_TOS_DATA\runs\batos_paper_v2\<crop_id>\` with a `run.json`.
