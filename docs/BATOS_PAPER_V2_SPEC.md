# BA-TOS paper-v2 re-implementation spec (preset `batos-paper-v2`)

Written 2026-11-04 (Phase 1a, documentation only; no code written).

**Authoritative paper text:** the camera-ready source
`C:\Users\USER\Downloads\2026_conf_sibgrapi_main_track_batos_felipe_dias_camera_ready\main.tex`
(2026-09-15). I checked its compiled `main.pdf` in the same folder against it (pp. 3–4). Text and
numbers are identical. Table numbering in the compiled PDF is I = frozen parameters,
II = crop-macro, III = pooled ablation, IV = parameter sensitivity, V = crop extent.
**There is no paper PDF or .tex in this repo's `paper/`.** `paper/figs/` holds only four PNGs
(fig2a, fig2b, fig3_crop02, fig4_crop05) from the earlier v4 draft. The camera-ready uses different figures (`fig1_pipeline_batos`,
`fig2_marker_diagnostic`, `fig5_visual_comparison_*`), so I could not compare against a PDF in `paper/`.

*Update 2026-11-04 (Phase 1b):* `paper/batos_sibgrapi_2026_camera_ready.pdf` has since been added. Its SHA-256
equals the camera-ready `main.pdf` (prefix `ed04bafd0bd72388`), so there is no difference from `main.tex`.
`docs/CLAUDE_CODE_BATOS_JOURNAL_PLAN.md` is **still absent**.

Companion documents: `docs/RECOVERED_SCRIPTS_INVENTORY.md` (script-level detail) and
`data/metadata/crops_sibgrapi.json` (crops). `docs/CLAUDE_CODE_BATOS_JOURNAL_PLAN.md` was **not
found** in the repo or anywhere on disk, so this spec does not use it.

---

## 1. Key finding: two pipelines are reported in one paper

The camera-ready mixes numbers from **two different pipelines**:

| | Legacy pipeline ("May-2026", README era) | E-series pipeline (E10 → E11 → E13 → E23) |
|---|---|---|
| Preprocessing | Gaussian flattening, σ = λ_xy = **15**, on the raw crop | GT-free scale estimate R → **opening-by-reconstruction**, radius max(8, round(1.5R)) → Gaussian flattening with **λ = 1** on that output |
| ToS candidates | `batos_tos_seed.py`, K (`--max-selected-nodes`) = **900** | `batos_tos_seed_true_dual.py`, K = **3600**, score weights listed in inventory §0 |
| Markers | `batos_background.py` + `batos_markers.py` (component markers, NMS 5/7, ≤ 420) | one point per internal component, priority score + 0.2·EDT + 0.05·MIN, NMS **6** |
| Foreground mask | Otsu **− 10** on robust uint8, closing ball(1), marker forcing | Otsu (no offset) on p1–p99 image, 6-conn closing ×1, fill holes, components ≥ 500, ∪ mask markers |
| Size filter | predictions < 500 removed | none in the watershed script; whether e09b applies one is unknown |
| Evaluation | `batos_instance_metrics.py`: Hungarian, IoU 0.50, border/area-ignore by overlap majority | `e09b_ignore_aware_audit.py`: **seed-ignore**, IoU 0.25/0.50/0.75 |
| Numbers in paper | **Table III** (pooled 0.980/0.967/0.974, FP 39, FN 77, ablations) | **Tables II, IV, V**, gate statistics, abstract crop-macro numbers |
| Table I entries that describe it | λ_xy = 15; "Otsu offset −10 / 1 voxel"; 500 | τs 0.35, K 3600, NMS 6, gate, IoU 0.25/0.50/0.75 |

### 1.1 The λ discrepancy: my reading, with evidence

**Reading:** the published crop-macro, sensitivity, crop-extent and gate numbers, which are all of
the Phase 1 fidelity targets except pooled 0.974, were produced with **opening-by-reconstruction
followed by λ = 1**, not λ_xy = 15. Table I's "λ_xy = 15" (and "Otsu offset −10") describe
the legacy pipeline. That pipeline produced only the pooled ablation in Table III. Table I is therefore a
hybrid of the two configurations.

Evidence:
1. **E13 code is explicit.** `e13_estimate_scale.py` hard-codes `"lambda_after_openrec": 1`.
   `e13_reconstruct_official_pipeline.py` calls the openrec script and then
   `run_batos_pipeline.sh <openrec> <out> 1 <min_size>`. `e13_validate_reconstruction.py` fails
   if a selector is `lambda15p0` and requires `flatten_lambda1p0`. The README says "Do not continue if one crop
   uses `lambda15p0`".
2. **E11, the frozen run before E13, already used openrec + λ1.** `run_e11_all_crops.sh` and
   `run_e11_holdout_suite.sh` read `01_preprocess_openrec_auto_minparam/*_flatten_openrec_r*_8u.tif`
   and `batos_openrec_lambda1_min500/01_preprocess/*_flatten_lambda1p0_8u.tif`. E13 only
   re-ran that preprocessing uniformly on all 12 crops. Every downstream number (E23c 0.823817 →
   0.906272; E24 sensitivity; E25 crop size) comes from these roots.
3. **Table III is the README result.** P 0.980, R 0.967, F1 0.974, FP 39, FN 77 appear verbatim in
   `README.md`, committed **2026-05-08**. That is months before E11 (the hotfix README names
   `e11_all_crops_frozen_20260725_*`). In May the pipeline was `scripts/run_second_scan_EFRGP01_01_selected_crops.sh`
   → `run_batos_pipeline.sh … 15 500`, i.e. λ = 15, plus `batos_instance_metrics.py --min-gt-area 1000
   --ignore-gt-border --thresholds 0.50`.
4. **The numbers cannot come from one pipeline.** Table III's FN = 2808 for "w/o boundary
   rejection" (TP = 0) and the strict-audit mean of 234.0 references × 12 crops = 2808 give the same
   valid-GT pool. BA-TOS_ToS has crop-macro mean 0.824 (sum 9.888). The ten non-crop05 crops lie in
   0.950–0.991 (Discussion), which sum to at least 9.50, so the two crop05 crops together have F1 ≤ 0.39 (one is 0.012).
   With roughly 230 valid GT per crop, those two crops alone would give hundreds of FN or FP in a pooled count. That is
   incompatible with FN 77 / FP 39. So the pooled 0.974 is **not** the pooled score of the
   BA-TOS_ToS route in Table II.
5. The paper's own Eq. (1) only says B_σ is estimated "at a spatial scale larger than the expected grain
   radius". Openrec with r = 1.5R fits that description as well as a Gaussian σ = 15 does.

Confidence: high that Tables II/IV/V and the gate used openrec + λ1. Medium on exactly how openrec
computes B, because `batos_preprocess_opening_reconstruction_v2.py` is missing.

**Consequence for Phase 1:** `batos-paper-v2` implements the E-series route (openrec + λ1). A second
preset, `batos-legacy-may2026` (λ = 15, legacy markers and mask, Hungarian evaluator), is the only
route that can target Table III / pooled 0.974. CLAUDE.md §2 lists pooled 0.974 as an expected
result of the same reproduction, so §2 and §4 need a note. **No frozen value is changed here.** This is
flagged for the author to decide (open question Q1).

---

## 2. Stage-by-stage spec (Algorithm 1 + §II of the paper)

Legend: **Exists** = recovered or in-repo code usable as reference or as-is. **Rewrite** = must be
written in `experiments/batos_v2/` (never edit `src/batos` or the recovered scripts).
"Fidelity" lists the checks from CLAUDE.md §4 Phase 1 plus stage-level checks I propose.

### 2.0 Crop extraction
- **Exists:** origins (`e13_build_12crop_manifest.py` CROP_MAP, `run_second_scan…sh`) are now in `data/metadata/crops_sibgrapi.json`.
- **Rewrite:** `batos_make_paired_crop.py` is missing. Read a 200³ window from the (Z,Y,X) uint16 BigTIFF with a memmap, without loading 6 GB. GT comes from `EFR_labels.7z` (extraction is a later, separate step).
- **Frozen:** size 200; origins as listed.
- **Fidelity:** axis order unverified. Compare our crop02 at (600, 800, 300) of EFRGP01_00 with the surviving `Desktop\batos codigos\lysandre_code\experiments\crops\crop_600_800_300_200.tif` (≈8.0 MB, i.e. 200³ uint8). Try both `[z:z+200, y:y+200, x:x+200]` and the x/z swap. Expect corr ≈ 1 for exactly one ordering after uint8 rescaling.

### 2.1 Conservative flattening (Alg. 1 line 1)
- **Exists:** λ step = `src/batos/batos_preprocess.py` (per-slice Gaussian σ = λ, I − B + median(B>0), p0.5–p99.5 → uint8). Scale estimate = `e13_estimate_scale.py` (complete, GT-free).
- **Rewrite:** opening-by-reconstruction flattening (the `_v2` script is missing). Proposed: B = reconstruction-by-dilation of the grayscale opening with ball(r) under I, with r = max(8, round(1.5R)); I_f = I − B + median(B), then robust uint8. Then apply `flatten_volume_u8(…, λ=1)` → **selector image**. Record the uint8 openrec image as the **watershed image**.
- **Frozen:** r = max(8, round(1.5R)); λ_after_openrec = 1; robust percentiles 0.5/99.5 (λ step). The paper's λ_xy = 15 applies to the legacy preset only.
- **Fidelity:** no recovered intermediate survives, so the only checks are downstream (internal mask non-empty; gate statistics; F1). Log R and r per crop and report them.

### 2.2 ToS MIN candidates, K = 3600 = 1800 + 1800, τs = 0.35 (line 2)
- **Exists (updated after the Q5 diff):** `origin/batos-ls-auto-tile-pilot:CodeReport_LysandreMacke/code/segmentation/batos/batos_tos_seed.py` (437 lines). It is the direct ancestor of the missing `batos_tos_seed_true_dual.py`: it has the **identical flag set** (height/area factors, dark percentile, support terms, 5 weights, min-score, max-selected-nodes, cleanup radii, component min/max), builds `hg.component_tree_tree_of_shapes_image3d`, classifies nodes MIN-/MAX-like, scores them, truncates to `max_selected_nodes`, and writes `*_tosconop_autodual_candidates.csv`. `e14_build_mintree_candidates.py` gives the recovered `candidate_diagnostics` conventions.
- **Port + one change:** port `batos_tos_seed.py` unchanged and add the "true dual" capacity split (best guess: rank MIN-like and MAX-like candidates separately and keep up to K/2 of each, instead of one global truncation to K). Use the Higra ToS on the 3D selector image (6-connectivity). Compute per-node area, parent–child contrast (height), darkness percentile ≤ 78, support terms (support area 20–45000, ≤ 16 steps up the ancestry) and the weighted score (w_height .22, w_delta .12, w_area .18, w_dark .18, w_support .30; area preference at log-position 0.42, width 0.55; area factor 22; area 2–2200; component 2–1500). Keep score ≥ τs. Rank, keep up to K/2 MIN-like and K/2 MAX-like nodes, and write the candidate CSV (`node, polarity, score, …`).
- **Frozen:** τs 0.35; K 3600 (1800/1800); area factor 22; area range 2–2200; all weights above.
- **Fidelity:** τs sweep on crop04 is flat at 0.985 / 0.978. K sweep 900/1800/3600/7200 gives F1@.50 0.510/0.833/0.985/0.998 (F1@.75 0.478/0.809/0.978/0.998). Retained ToS markers per crop fall in **53–405**. MIN and DUAL marker sets must be **identical**.

### 2.3 Spatialization and exclusive supports (Alg. 1 line 5, input to the gate)
- **Exists:** output schema only (`candidate_diagnostics_no_gt.csv`; inventory §0).
- **Rewrite:** `e10_spatialize_candidates_no_gt.py`. Per selected node: reconstructed region; representative point (`point_z/y/x`); `has_reconstructed_region`; `region_touches_volume_face`. `region_voxels` = |R_i|, where each leaf voxel goes to its **nearest selected ancestor** (paper Eq. 5).
- **Fidelity:** both crop05s give m_v = 2, q_v = 3; all other 10 crops trigger nothing, with the closest at m_v = 17, q_v = 43. The median region over crop04 ≈ 43 voxels (E23B README).

### 2.4 Boundary rejection, C ∩ ∂Ω ≠ ∅ (line 3)
- **Exists:** `src/batos/batos_background.py` → `internal_dark_candidates.tif`. Dark threshold = min(Otsu − 10, P40, 135) on the selector image. Border-seeded 6-conn propagation; internal = dark ∧ ¬external. In the E-series, a ToS candidate survives only if its point lies in a 6-connected `internal_dark` component (`e11_build_markers_ablation.py`). So the operational rule is "**the candidate's dark connected region does not reach any face**". That is equivalent to C ∩ ∂Ω = ∅ applied to the dark level-set component, not the ToS node itself.
- **Rewrite:** nothing, provided `batos_background.py` in the repo equals the Linux version (verify by reading; see Q5). It needs `--seed`, but the seed does not affect `internal_dark`. Pass an empty seed or the new raw seed, and record which.
- **Frozen:** `--dark-percentile 40 --otsu-offset -10 --dark-threshold-max 135 --dark-connectivity 1` (from `run_batos_pipeline.sh`).
- **Fidelity:** ablation "w/o boundary rejection" → F1 0.000 (Table III, legacy). In v2, check that the internal mask is non-empty on all 12 crops (E13 validation rule).

### 2.5 Internal markers + NMS 6 (line 4)
- **Exists:** `e11_build_markers_ablation.py` is complete and deterministic. Port it as-is: one point per internal component, priority = minmax(score) + 0.20·log1p(EDT_max(comp))/log1p(max) + 0.05·[MIN], greedy NMS with r = 6 (strict `<` r²), single-voxel markers in uint16.
- **Also needed:** `baseline_existing` markers (`src/batos/batos_markers.py`, legacy flags `--area-min 3 --area-max 900 --max-dim 30 --max-elongation 8 --min-fill 0.005 --center-dark-max 205 --min-ring-contrast -2 --nms-min-dist 7 --max-markers 420`). They are **not** a compared method in v2, but they **enlarge the frozen watershed mask** (`--mask-marker`).
- **Frozen:** NMS 6; EDT weight 0.20; MIN bonus 0.05; polarity MIN; connectivity 6.
- **Fidelity:** NMS 4/6/8/10 on crop04 gives F1@.50 0.976/0.985/0.985/0.985.

### 2.6 Support-collapse gate: median ≤ 3 ∧ q75 ≤ 4 → EDT-Otsu (lines 5–6, 8)
- **Exists:** `e23b_degeneracy_gate_audit.py` and `e23c_build_selective_marker_route.py` are complete. The gate takes `region_voxels` of the **selected** ToS MIN marker nodes (merge on `node`), median and q75 via pandas `quantile` (linear). EDT-Otsu fallback = `e11_build_control_markers.py --mode otsu_edt --nms-radius 6 --match-count-from <ToS markers>` on the **watershed (openrec) image**. Its count equals the ToS marker count of the same crop.
- **Rewrite:** a port only.
- **Frozen:** 3 / 4; fallback EDT-Otsu with matched count.
- **Fidelity:** exactly {EFRGP01_00_crop05, EFRGP01_01_crop05} trigger with (2, 3). The 10 non-gated segmentations must be bitwise identical between the ToS and sel routes.

### 2.7 Mask: "Otsu −10 + closing 1" (line 7)
- **Exists:** `e11_run_fixed_mask_watershed.py` (complete).
- **Discrepancy (second, after λ):** the recovered E-series mask is **plain Otsu** (`>` thr on the p1–p99 normalized openrec image), with **no −10 offset**, then 6-conn closing ×1, fill holes, components ≥ 500, OR the ToS MIN and baseline marker voxels. The "Otsu offset −10 / closing 1 voxel" in Table I matches the legacy `batos_watershed_distance.py --threshold-offset -10 --close-radius 1`. In the E-series, "−10" appears only in the dark threshold of the boundary step (2.4). v2 follows the E-series code, because the fidelity targets come from it.
- **Fidelity:** mask SHA-256 is logged per crop, and the gate must not change the mask.

### 2.8 Watershed on −EDT (line 9)
- **Exists:** same script. Relief = −(EDT normalized to p1–p99 within the mask); `skimage.segmentation.watershed(connectivity=1, mask, watershed_line=False)`; markers clipped to the mask and relabeled with 6-conn.
- **Note:** the paper credits Higra for the watershed primitives, but the recovered code uses scikit-image. v2 uses scikit-image.

### 2.9 Size filter 500 (line 10)
- **Exists:** nothing in the E-series watershed. The legacy `filter_predicted_labels_by_size.py` is missing (trivial). E13 also derives a per-crop `min_predicted_size = max(500, round(0.1·V̂))`, but passes it only to the legacy runner.
- **Rewrite:** `remove predicted labels < 500 voxels` and relabel.
- **Open:** whether the published v2 numbers applied it (in e09b?) and with 500 or the per-crop value (Q3). Run v2 **with and without** it, and report both against the targets.

### 2.10 Seed-ignore evaluation at IoU 0.25 / 0.50 / 0.75
- **Exists:** `e12_boundary_margin_audit.py` (strict protocol with max-cardinality bipartite matching) is the nearest recovered evaluator. `src/batos/batos_instance_metrics.py` is legacy.
- **Rewrite:** `e09b_ignore_aware_audit.py`. Valid GT = ≥ 1000 voxels and not touching any face (bbox start ≤ 0 or stop ≥ size). A prediction is **ignored** iff its seed marker lies inside an excluded GT label. A prediction seeded in background is evaluated (it is an FP if unmatched). Matching: maximum bipartite matching on IoU ≥ thr (as in e12). Output columns `{tp,fp,fn,precision,recall,f1}_iou{025,050,075}` with `evaluation_protocol ∈ {seed_ignore, strict}`.
- **Open:** how predictions that overlap ignored GT, but whose seed is elsewhere, contribute to IoU (Q4). Verify the evaluator on `sample_data/` (self-match → F1 = 1).
- **Fidelity:** dev crop04 F1@0.50 ≈ 0.985 (F1@.75 0.978). Crop-macro mean F1@0.50: BA-TOS_ToS 0.824 ± 0.367 (median 0.985, min 0.012); BA-TOS_sel 0.906 ± 0.210 (min 0.260). Mean F1@0.75 0.817 → 0.899. Ten gate-negative crops in [0.950, 0.991]. Crop-extent check: 128/160/200 → 0.851/0.967/0.985.

### Baselines (for later; not Phase 1)
Table II's "Otsu+WS" (0.830 macro, min 0.260) shares its minimum with BA-TOS_sel (0.260). That suggests Otsu+WS = `edt_otsu` under the same frozen mask (it is identical to the fallback on crop05). "Min-tree-core+WS" (0.768) may be `mintree_e14` or the Macke baseline (`run_macke_mintree_one_crop.sh`, missing). Q6.

---

## 3. Fidelity acceptance table (Phase 1b/1c)

| Check | Target | Source |
|---|---|---|
| crop04-dev F1@0.50 / @0.75 | 0.985 / 0.978 | Table IV |
| Gate fires | only both crop05, (m_v, q_v) = (2, 3); closest other (17, 43) | §IV-B |
| K sweep, crop04 | 0.510 / 0.833 / 0.985 / 0.998 | Table IV |
| τs sweep, crop04 | flat 0.985 | Table IV |
| NMS sweep, crop04 | 0.976 / 0.985 / 0.985 / 0.985 | Table IV |
| Crop-macro F1@0.50 | ToS 0.824 ± 0.367, sel 0.906 ± 0.210 | Table II |
| Crop-macro F1@0.75 | 0.817 → 0.899 | Table II |
| Retained ToS markers | 53–405 per crop | §II-B |
| Crop extent | 0.851 / 0.967 / 0.985 | Table V |
| Pooled core 0.974 (P .980 R .967 FP 39 FN 77) | **only via `batos-legacy-may2026`** | Table III, README |

Any deviation is reported, not tuned away (CLAUDE.md rule 4).

---

## 3b. Author decisions (2026-11-04) and best-guess register

Decisions: v2 = E-series pipeline. The legacy λ15 route becomes the documented preset `batos-legacy-may2026`. Pooled F1 is recomputed with v2 and not forced to 0.974. Table II baselines: Min-tree-core+WS = `mintree_e14`; Otsu+WS = `edt_otsu` matched-count. Seed-ignore follows the paper's §III-C literally. The prediction size filter is a switch, and Phase 1c reports both settings.

**Best-guess register.** Each item is reconstructed, not recovered. Each must be YAML-configurable and logged per run.

| ID | Step | Best guess | Basis | Config key |
|---|---|---|---|---|
| BG1 | Scale estimate R̂ | Port `e13_estimate_scale.py` exactly (p1–p99 uint8, Otsu ≥, closing ball(1), fill holes, 5³ EDT peaks ≥ max(2, P75), median of top 1500, fallback P95, floor 3) | recovered code | `scale.method: e13` |
| BG2 | Openrec radius | r = max(8, round(1.5·R̂)), per crop | recovered code | `openrec.radius: auto` (or an int), `openrec.radius_factor: 1.5`, `openrec.radius_min: 8` |
| BG3 | Openrec operator | grayscale **opening by reconstruction**: erosion by ball(r), then reconstruction-by-dilation under I (skimage `morphology.reconstruction`) | script name + E13 docstring "opening by reconstruction" | `openrec.structuring_element: ball`, `openrec.connectivity: 6` (6 vs 26 unknown) |
| BG4 | Openrec background use | The output is named `flatten_openrec`, so the reconstruction is treated as a **background estimate**: I_f = I − B + median(B) (paper Eq. 1), clipped ≥ 0 → robust uint8 (p0.5/p99.5 on >0 voxels, as in `batos_preprocess.py`). Alternative to test if fidelity fails: I_f = B (the openrec image itself, i.e. a smoothing rather than a flattening). | paper Eq. 1, naming | `openrec.mode: subtract_background \| reconstructed_image` |
| BG5 | Polarity of the opening | Opening (removes bright structures smaller than r). The grains are bright, so the opening of a grain-sized structure is the dark background. Closing variant kept as an option. | grains bright in E11 mask (`image > otsu`) | `openrec.operation: opening \| closing` |
| BG6 | Selector flattening | `batos_preprocess.flatten_volume_u8(openrec_u8, lambda_xy=1, sigma_z=0)` | E13 runner + repo code | `selector.lambda_xy: 1.0` |
| BG7 | Input dtype for openrec | Raw uint16 crop (E13 manifest passes the `*_gray_crop.tif` window from `batos_make_paired_crop.py`, which writes the raw dtype). Robust-u8 conversion happens after. | pilot-branch crop script | `openrec.input: raw` |
| BG8 | True-dual split | K/2 MIN + K/2 MAX by separate ranking | paper §II-B text | `tos.k_total: 3600`, `tos.k_split: half` |
| BG9 | Exclusive support | Each leaf voxel goes to its nearest selected ancestor; `region_voxels` = count | paper Eq. 5 | — |
| BG10 | Size filter | `size_filter.enabled: true/false`, `size_filter.min_voxels: 500` (per-crop `max(500, 0.1·V̂)` as a third option, not run unless asked) | E13 params, paper | `size_filter.*` |
| BG11 | Seed of a prediction (seed-ignore) | The marker voxel(s) of the watershed basin. Ignored iff **all** seed voxels fall in excluded GT labels; evaluated if any seed voxel is in background or valid GT. Matching = maximum bipartite on IoU ≥ thr over the domain minus excluded-GT voxels (as in e12). | paper §III-C + e12 | `eval.protocol: seed_ignore` |

Draft YAML layout for `configs/batos_paper_v2.yaml` (to be written in Phase 1c, not now):
```yaml
preset: batos-paper-v2
scale: {method: e13}
openrec: {radius: auto, radius_factor: 1.5, radius_min: 8, operation: opening, mode: subtract_background, connectivity: 6, input: raw}
selector: {lambda_xy: 1.0, sigma_z: 0.0}
boundary: {dark_percentile: 40, otsu_offset: -10, dark_threshold_max: 135, connectivity: 1}
tos: {k_total: 3600, k_split: half, min_score: 0.35, area_factor: 22.0, area_min: 2, area_max: 2200}
markers: {polarity: min, nms_radius: 6, edt_priority_weight: 0.20, min_priority_bonus: 0.05}
gate: {median_max: 3, q75_max: 4, fallback: edt_otsu_matched_count}
mask: {otsu_offset: 0, closing_iterations: 1, fill_holes: true, min_component: 500, union_markers: [tos_min, baseline_existing]}
watershed: {surface: neg_edt_p1p99, connectivity: 1}
size_filter: {enabled: true, min_voxels: 500}
eval: {protocol: seed_ignore, iou: [0.25, 0.50, 0.75], min_gt_voxels: 1000, exclude_face_touching: true}
```

## 3c. Controlled reconstruction experiments on EFRGP01_00_crop04 (2026-11-04) and frozen choice

Authorised scope: development crop only. Only reconstructed steps may change. K, τs, NMS, the boundary rule, the gate,
the mask and the watershed are frozen (enforced by `scripts/run_v2_control.py`). Target: F1@0.50 ≈ 0.985, F1@0.75 ≈ 0.978.
All runs are under `C:\BA_TOS_DATA\runs\batos_paper_v2_controls\<name>\result.json`. All F1 values are seed-ignore,
and identical with the size filter at 500 and off, except H4 (500 shown).

| Run | Factor changed (one at a time vs v2 defaults) | MIN pool | Markers | Valid grains seeded /272 | F1@0.50 | F1@0.75 |
|---|---|---|---|---|---|---|
| v2 main / C2a | none (openrec r15 subtract → λ1 selector) | 1,938 | 255 | 172 | 0.559 | 0.491 |
| C1 | original `batos_tos_seed.py` (pilot branch, unmodified) on the same image | 1,938 | — | — | — | — |
| C2b | selector = main `batos_preprocess.py` λ15 on raw (May image) | 858 | 1 | 0 | 0.000 | 0.000 |
| C2c | selector = λ1 on raw crop, no openrec | 1,945 | 261 | 178 | 0.609 | 0.538 |
| C3 | openrec 8-bit stretch p1–p99 | 1,902 | 250 | 167 | 0.542 | 0.501 |
| H1 | openrec mode = top-hat, max(I − B, 0) | 1,906 | 256 | 176 | 0.571 | 0.527 |
| H2 | 3D Gaussian σ = 1 before the opening | 1,632 | 111 | 65 | 0.137 | 0.125 |
| H3 | opening on the robust-u8 crop | 1,885 | 261 | 173 | 0.575 | 0.512 |
| H4 | 2D per-slice opening, disk r | 2,421 | 261 | 168 | 0.532 | 0.431 |
| H5 | division I/B·median(B) | 1,938 | 255 | 172 | 0.559 | 0.491 |

Findings:
1. **The port is exact (C1).** The original selector and `batos_v2.tos_select` give the identical 1,938-node MIN list, same order, zero difference in every score component.
2. **λ15 is incompatible with the E-series boundary step (C2b).** Almost all dark structure is border-connected, one marker survives, and the gate fires. This independently confirms why E13 rejected `lambda15p0`.
3. **The opening-by-reconstruction step is not the cause of the gap.** At r ≥ the grain radius the background B is nearly flat, so subtraction, division and top-hat are near-identity (C2a ≈ C2c ≈ H1 ≈ H5). Every openrec variant leaves the MIN pool at about 1.6k–2.4k and seeded grains at ≤ 178/272. The paper's K sweep (0.985 → 0.998 from K = 3600 to 7200) implies an original MIN pool well above 3,600.
4. **Frozen choice:** the a-priori BG defaults are kept (subtract_background, raw, 3D, no pre-smoothing, 0.5/99.5). Picking the best of several near-equal dev-crop variants (C2c, H3, H1: +0.016 to +0.05) would be tuning without a supported mechanism.
5. **Fidelity is not reached: v2 gives F1@0.50 = 0.559 on crop04 vs 0.985.** The remaining suspects are outside the openrec step and were not authorised in this round:
   - (a) The true-dual selector's MIN scoring/capacity, e.g. K applied to raw nodes before τs, or MIN statistics such as `avg_min_area` and `max_min_height` computed differently.
     On this image `max_min_height` = 10, so `h_score` saturates at 1 for every node and τs = 0.35 is decided by the delta, area, darkness and support terms.
   - (b) The representative-point rule (BG12): 62 % of MIN points fall on inter-grain background and fail internal eligibility.
   - (c) The lost e10 spatialization's eligibility, which may project candidates into internal components rather than requiring the point inside one (E11 carries a `projection_distance` column, always 0 in the recovered builder).

## 4. Environment plan

**Status 2026-11-04: implemented (Phase 1b).** The venv `C:\BA_TOS_DATA\envs\batos-py312` (Python 3.12.10) is created and frozen to `environments/batos-py312.txt`. Higra 0.6.13 has a cp312 wheel, and a 3D ToS smoke test passed. The plan as originally written follows.


Current machine (2026-11-04): Python **3.14.4** (py launcher default) and **3.12.10**
(`C:\Users\USER\AppData\Local\Programs\Python\Python312\python.exe`). Both have numpy, scipy and
CPU-only torch. Neither has tifffile, py7zr, higra, scikit-image, pandas or pyyaml. The RTX 5060
(8 GB, driver 617.14) is not visible to the installed torch builds.

Proposal:
- **Python 3.12.10**. Higra's published wheels are the binding constraint, and 3.14 wheels are unlikely. Fall back to 3.11 only if `pip install higra` fails on 3.12.
- **venv:** `C:\BA_TOS_DATA\envs\batos-py312` (outside the repo, so it is never committed).
- **Packages** (pin exact versions in `configs/env_batos_py312.txt` after the first install and record them in the log): `numpy scipy scikit-image tifffile py7zr higra pandas pyyaml` (+ `pillow matplotlib` from `requirements.txt` for figures).
- No Cellpose and no CUDA torch in this env (Phase 2 gets its own env).
- Run all data-reading Python with `-I`. Write outputs under `C:\BA_TOS_DATA\` (data) or `results/` (git-ignored).

Proposed commands, **for Phase 1b, not run**:
```
"C:\Users\USER\AppData\Local\Programs\Python\Python312\python.exe" -m venv C:\BA_TOS_DATA\envs\batos-py312
C:\BA_TOS_DATA\envs\batos-py312\Scripts\python.exe -m pip install --upgrade pip
C:\BA_TOS_DATA\envs\batos-py312\Scripts\python.exe -m pip install numpy scipy scikit-image tifffile py7zr higra pandas pyyaml pillow matplotlib
C:\BA_TOS_DATA\envs\batos-py312\Scripts\python.exe -m pip freeze > configs\env_batos_py312.txt
```
