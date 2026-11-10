# BA-TOS journal extension — experiment log

Append-only (CLAUDE.md rule 8). Newest entries at the bottom.

---

## 2026-11-04 — Session 1: Phase 0 + Phase 1a (read-only inspection, documentation)

Machine: native Windows 11 Pro 10.0.26200 (PowerShell + Git Bash). Repo `C:\Users\USER\batos-3d-grain-segmentation`,
branch `main`, HEAD `478a4fc`. No code was run on scientific data. Nothing was installed. The 7z was not extracted.
This session was interrupted once by a connection error and resumed. CLAUDE.md was then copied into
the repo root (byte-identical to `Downloads\CLAUDE.md`), and the patch zip was extracted by the user.

### Read
- `CLAUDE.md` (repo root, identical to the Downloads copy).
- **`docs/CLAUDE_CODE_BATOS_JOURNAL_PLAN.md`: NOT FOUND.** It is not in the repo, not in Downloads, Desktop or Documents, and not under `C:\BA_TOS_DATA`. Phase 0/1a went ahead without it.
- Repo: `README.md`, `requirements.txt`, `.gitignore`, every file in `src/batos/` and `scripts/`, `sample_data/` metadata.
- All 94 files in `C:\BA_TOS_DATA\linux_import\linux_patches\`. I verified them against `Downloads\batos_linux_patches_v2.zip` by SHA-256: 0 mismatches. In the first half of the session (before extraction) I read them in memory from the zip.
- Camera-ready source `Downloads\2026_conf_sibgrapi_main_track_batos_felipe_dias_camera_ready\main.tex` (2026-09-15) and its `main.pdf` pp. 3–4, which agree.
- `paper/` holds only `figs/` (4 PNGs from the v4 draft). It has no paper PDF/.tex.

### CLAUDE.md §1 files: existence and sizes
| Item | Status | Size |
|---|---|---|
| `src/batos/batos_instance_metrics.py` | present | 7,009 B |
| `src/batos/count_tiff_labels.py` | present | 828 B |
| other `src/batos/*.py` (background, evaluate_against_gt, markers, preprocess, watershed_distance) | present | 16,569 / 11,126 / 11,825 / 2,827 / 5,577 B |
| `scripts/` (run_batos_pipeline.sh, run_second_scan_EFRGP01_01_selected_crops.sh, make_two_scan_final_comparison.py) | present | 6,453 / 4,661 / 5,976 B |
| `sample_data/` (2 synthetic 96³ crops + labels, previews, csv, md) | present | 5.4 MB total |
| `paper/figs/` | present, 4 PNGs | — |
| `requirements.txt` | present (numpy scipy scikit-image tifffile pillow matplotlib pandas py7zr higra) | 80 B |
| `~/Documents/BA_TOS_DATA/EFRGP01_0{0,1}.tif` (CLAUDE.md path) | **absent**; see the volumes below | — |

### Raw volumes (TIFF headers parsed with a stdlib reader; no pixel loading)
| Copy | Path | Bytes | Format |
|---|---|---|---|
| **canonical** 00 | `C:\Users\USER\Downloads\PhDUpload\PhDUpload\EFR\Greyscale\EFRGP01_00.tif` | 6,272,409,728 | BigTIFF, uint16, 1600 pages of 1400×1400 → (Z,Y,X) = (1600,1400,1400), uncompressed |
| **canonical** 01 | `C:\Users\USER\Downloads\PhDUpload\PhDUpload\EFR\Greyscale\EFRGP01_01.tif` | 6,272,409,728 | same |
| secondary 00 | `C:\Users\USER\Desktop\BATOS_DATA\EFRGP01_00.tif` | 6,272,409,776 (+48 B) | same shape and dtype |
| secondary 01 | `C:\Users\USER\Desktop\BATOS_DATA\EFRGP01_01.tif` | 6,272,409,728 | same shape and dtype; **differs from the Downloads copy starting at byte 385** (header region) |

Decision: **the Downloads copies are canonical.** Pixel equality between the copies is **not verified**. It is a candidate check for Phase 1b (compare a few z-planes).
GT: `Downloads\PhDUpload\PhDUpload\EFR\Label\EFR_labels.7z` (2,623,879,422 B; 55 volumes EFRGP01_00..54 per `Desktop\batos codigos\lista_EFR_labels_py7zr.txt`). The sibling folders `extracted_one/` and `extracted_EFRGP01_01/` are **empty**.

### Hardware / Python
- CPU AMD Ryzen 7 5700X (8C/16T); RAM 31.9 GB; GPU NVIDIA GeForce RTX 5060, 8151 MiB, driver 617.14.
- C: has 115.3 GB free.
- Python 3.14.4 (`AppData\Local\Python\pythoncore-3.14-64`, py default): numpy 2.4.6, scipy 1.17.1, torch 2.12.0+cpu (CUDA False). Missing: tifffile, py7zr, higra, scikit-image, pandas, pyyaml.
- Python 3.12.10 (`AppData\Local\Programs\Python\Python312`): numpy 2.5.0, scipy 1.18.0, torch 2.12.1+cpu (CUDA False). Missing: the same six.
- Conclusion: no current interpreter can run any BA-TOS script. The environment plan is in `docs/BATOS_PAPER_V2_SPEC.md` §4 (Python 3.12, venv `C:\BA_TOS_DATA\envs\batos-py312`).

### Findings
1. **Two pipelines in one paper** (detail in the spec §1). Tables II/IV/V and the gate come from the E-series route (opening-by-reconstruction r = max(8, round(1.5R)) → Gaussian λ = 1 → ToS true-dual K3600 → point markers NMS 6 → fixed Otsu mask (no offset) → −EDT watershed → seed-ignore evaluation). Table III (pooled 0.974, FP 39, FN 77) equals the README of 2026-05-08, which comes from the legacy λ = 15 pipeline. Arithmetic check: 10 gate-negative crops in [0.950, 0.991] and mean 0.824 force the two crop05 crops to have summed F1 ≤ 0.39. That is incompatible with pooled FN 77 / FP 39 over 2808 valid GT. **Table I's λ_xy = 15 and "Otsu offset −10 / closing 1" describe the legacy route, not the one behind the Phase 1 fidelity targets.**
2. The boundary rule is implemented operationally as "the ToS point lies in a 6-connected `internal_dark` component", where `internal_dark` = dark voxels (≤ min(Otsu − 10, P40, 135)) not border-connected. It is not a test on the ToS node itself.
3. The EDT-Otsu fallback uses a **matched marker count** (equal to the crop's ToS count), with NMS 6, on the openrec image.
4. Missing core code: `batos_tos_seed_true_dual.py`, `e10_spatialize_candidates_no_gt.py`, `e09b_ignore_aware_audit.py`, `batos_preprocess_opening_reconstruction_v2.py`, `audit_true_batos_e07.py`, the Linux `run_batos_pipeline.sh` + `batos_tos_seed.py` + `filter_predicted_labels_by_size.py`, and `batos_make_paired_crop.py`. Fully recovered and portable: E11 marker builder, control markers, fixed-mask watershed, E12 strict evaluator, E13 scale estimator, E23b/c gate and router.
5. The repo's `scripts/run_batos_pipeline.sh` line 33 has a literal `\n` inside `echo` (cosmetic only; not changed).
6. `data/` is in `.gitignore`, so `data/metadata/crops_sibgrapi.json` will not be committed unless that rule is narrowed. CLAUDE.md §2b is dated "found 2026-11-05", one day after today's date (2026-11-04).

### Files written (docs/metadata only)
- `data/metadata/crops_sibgrapi.json`: 12 crops, origins, roles, `gate_expected`, `gt_path = C:\BA_TOS_DATA\gt\EFRGP01_0X_cropNN_labels.tif` (future), canonical `gray_path` = the Downloads copies, with Desktop as `gray_path_alt`.
- `docs/RECOVERED_SCRIPTS_INVENTORY.md`
- `docs/BATOS_PAPER_V2_SPEC.md`
- `docs/EXPERIMENT_LOG.md` (this file)

No existing code, results or frozen artefacts were modified.

---

## 2026-11-04 — Session 2: author decisions + Phase 1b (environment and data)

### Inputs / status
- Author decisions received (v2 = E-series; legacy preset for λ15/pooled 0.974; Q2/Q3 as best-guess or switch; Q4 = paper §III-C; Q5 diff; Q6 baselines = mintree_e14 and edt_otsu; Q7 via the Lysandre crop; Q8 via pixel SHA-256).
- **The instruction message was truncated** at "4. Run steps 2-3 only after Q7 is settled. Sanity check: the". The sanity check text is unknown. I assumed it is the valid-GT count against the paper (total 2808).
- `paper/batos_sibgrapi_2026_camera_ready.pdf` is present and byte-identical to the camera-ready `main.pdf` (SHA-256 prefix `ed04bafd0bd72388`).
- **`docs/CLAUDE_CODE_BATOS_JOURNAL_PLAN.md` is still absent** although it was reported as in place. It was not read.

### Docs updated
- `CLAUDE.md` §2 (flattening and mask split into v2 / legacy rows, size-filter switch, seed-ignore rule, v2 targets vs legacy pooled 0.974, baselines, a note on Table I), §2b (resolution note, pilot-branch recovery, canonical data, Desktop copy not to be used), §4 Phase 1 row (E-series, sub-phases 1a/1b/1c, both size-filter settings, pooled recomputed).
- `docs/BATOS_PAPER_V2_SPEC.md`: selector is now "port `batos_tos_seed.py` + true-dual split"; added §3b (decisions, best-guess register BG1–BG11, draft YAML); §4 marked implemented.
- `docs/RECOVERED_SCRIPTS_INVENTORY.md` §8: Q5 diff, the pilot-branch file list, and historical GT counts.
- `data/metadata/crops_sibgrapi.json`: canonical gray paths, pixel SHA-256 of all four copies, `gt_full_path`, `crop_convention` (axis_order `zyx`, **slicing_status `unverified`**, evidence), per-crop `gt_path`/`gray_crop_path` under `C:\BA_TOS_DATA\crops\<crop_id>\`, and `historical_gt_valid_legacy`.

### Q5: src/batos vs origin/batos-ls-auto-tile-pilot
Four modules are identical and three differ only by docstrings, shebangs or imports removed in the May "Clean up" commits. **`main:src/batos/batos_evaluate_against_gt.py` is broken (no imports)**; left untouched. The branch also contains `batos_tos_seed.py` (the ancestor of the missing true-dual selector, same flags), `batos_make_paired_crop.py`, `filter_predicted_labels_by_size.py`, `baseline_otsu_distance_watershed.py`, and historical result tables.

### Q8: pixel SHA-256 (uint16 little-endian, C order, streamed by memmap in 50-plane chunks, ~10 s/volume)
| Copy | pixel SHA-256 |
|---|---|
| Downloads EFRGP01_00 (canonical) | `4c6bfbf26f38ed37b51ba13e7d31076f9426397b74ec60a72a9ca0ef89ca0d26` |
| Desktop EFRGP01_00 | `72b4e1a2dc132b7a81bb206af6adf396df2b1f95a5d0da94162bb5091d75a612` |
| Downloads EFRGP01_01 (canonical) | `a20509ac2f4699550e175c89de6dfe97abb3f1476b6f2708bee265fe2e46bce2` |
| Desktop EFRGP01_01 | `8ee4d23f25a03a99e207ef7ef973e637214b7b296c22b0d828b3513e9fde4d72` |

**The pixel data differs for both scans.** At z = 300/800/1300 every voxel differs. The copies are aligned (corr 0.81/0.93/0.93 at equal z; 0.915 → 0.89 for z ± 1, 2) but the Desktop intensities are rescaled (mean ≈ 2.5k vs 13k) and planes z = 0 and z = 1599 are all zero. So the Desktop copy is a differently processed version. Downloads stays canonical; the historical Linux runner also read `/mnt/windows_p3/Users/USER/Downloads/PhDUpload/...`. A truncated third copy, `Downloads\Datas_contacts_points_test-001\EFRGP01_00.tif` (626 MB, unreadable header), was ignored.

### Phase 1b step 1: environment (DONE)
- `C:\BA_TOS_DATA\envs\batos-py312` created from Python 3.12.10. Installed: numpy 2.5.4, scipy 1.18.1, scikit-image 0.26.0, tifffile 2026.10.15, py7zr 1.1.4, **higra 0.6.13 (cp312 wheel available)**, pandas 3.1.1, PyYAML 6.0.3, plus dependencies.
- Frozen to `environments/batos-py312.txt`.
- Smoke test: `hg.component_tree_tree_of_shapes_image3d` on an 8×9×10 uint8 volume built 1149 nodes with 720 leaves (= voxel count) on Python 3.12.10.

### Q7: axis order (NOT settled)
- Script: scratchpad `q7_axis.py` (run with `-I`). For each of the 6 orderings of origin (600, 800, 300) onto the volume axes, it cut a 200³ block and compared it with the Lysandre crop under 6 transposes, raw and robust-u8.
- Result: **no match** on the Downloads or Desktop EFRGP01_00, or on EFRGP01_01. The best |corr| was 0.016, the same as unrelated texture. `crop_600_800_300_200.tif` (200³ uint8, tifffile-written, range 30–255) is not a crop of these scans at that origin. This method cannot settle Q7.
- Code evidence favors `zyx`: the original `batos_make_paired_crop.py` uses `arr[z0:z1, y0:y1, x0:x1]`.
- Proposed decisive check: the historical per-crop `gt_valid` for all 12 crops (sum **2808** = the paper's pool) is now recorded. After extracting GT, `make_crops.py --probe-axis-orders` reports the valid-GT count under all 6 orderings without writing anything. The correct ordering should reproduce 265/216/272/129/271/263 and 262/213/264/128/265/260.

### Phase 1b steps 2–3: scripts written, NOT run (blocked on Q7 per instruction)
- `scripts/extract_gt_labels.py`: py7zr extracts only `EFRGP01_00.tif` and `EFRGP01_01.tif` into `C:\BA_TOS_DATA\gt_full\`, reports shape, dtype and axes from the header, and writes `extract_report.json`. It refuses to overwrite without `--force`.
- `scripts/make_crops.py`: memmap only. Writes `C:\BA_TOS_DATA\crops\<crop_id>\{gray.tif (raw uint16), labels.tif}` and `data/metadata/crops_manifest.csv` (pixel and file SHA-256, shape, dtype, min/max, n_gt_labels, n_gt_touching_face, n_gt_small_lt1000, n_gt_valid). It **refuses to run unless `crop_convention.slicing_status == "verified"`**. `--probe-axis-orders` writes nothing.
- Both pass `py_compile`. Neither has been executed on data.

---

## 2026-11-04 — Session 3: plan read, GT extracted, Q7 settled, 12 crops cut

### Read
- `docs/CLAUDE_CODE_BATOS_JOURNAL_PLAN.md` (12,046 B), read in full. It is consistent with CLAUDE.md. Its "First Action" deliverable, `docs/BA_TOS_BASELINE_REPRODUCTION.md`, **does not exist yet**. The recovered-scripts inventory and the paper-v2 spec cover most of its content, but the dedicated document is still owed.
- Author answers: the sanity check is valid GT Σ = 2808, with the per-crop historical counts as the acceptance test. GT extraction and the axis probe were approved.

### GT extraction (`scripts/extract_gt_labels.py`, venv batos-py312, `-I`)
- `EFRGP01_00.tif` and `EFRGP01_01.tif` were extracted from `Downloads\PhDUpload\PhDUpload\EFR\Label\EFR_labels.7z` into `C:\BA_TOS_DATA\gt_full\` (35 s each). Both are uint16, (1600, 1400, 1400), tifffile axes `QYX` (= Z,Y,X), 6,272,409,776 B and 6,272,409,728 B. Report: `C:\BA_TOS_DATA\gt_full\extract_report.json`.
- Example: z = 800 of 00 has 616 distinct label values, max 22,321, 87 % zeros.

### Q8 correction (supersedes the Session 2 interpretation)
`cmp` shows `Desktop\BATOS_DATA\EFRGP01_00.tif` and `EFRGP01_01.tif` are **byte-identical to the extracted GT label volumes**. Their pixel SHA-256 (`72b4e1a2…`, `8ee4d23f…`) equal the GT hashes. They are not a differently processed greyscale copy. The rescaled "intensities", zeroed end planes and z-increasing means seen in Session 2 are label IDs. Downloads PhDUpload greyscale remains the only greyscale source. CLAUDE.md §2b and `crops_sibgrapi.json` are corrected.

### Q7: axis order SETTLED = `zyx`
`scripts/make_crops.py --probe-axis-orders` (writes nothing), as labels/valid per crop:

| crop | zyx | zxy | yzx | yxz | xzy | xyz | historical |
|---|---|---|---|---|---|---|---|
| 00_crop02 | 650/**265** | 641/260 | 418/182 | 428/193 | 229/88 | 247/96 | 265 |
| 00_crop03 | 473/**216** | 477/219 | 76/12 | 87/17 | 416/192 | 434/188 | 216 |
| 00_crop04 | 649/**272** | 649/272 | 663/270 | 639/267 | 663/270 | 639/267 | 272 |
| 00_crop05 | 304/**129** | 300/127 | 658/267 | 653/263 | 431/195 | 446/199 | 129 |
| 00_crop06 | 561/**271** | 542/265 | 510/239 | 501/235 | 644/266 | 654/270 | 271 |
| 00_crop07 | 627/**263** | 637/259 | 645/265 | 658/271 | 621/270 | 620/268 | 263 |
| 01_crop02 | 659/**262** | 651/261 | 422/182 | 440/190 | 228/87 | 252/97 | 262 |
| 01_crop03 | 470/**213** | 480/218 | 75/12 | 92/18 | 417/190 | 440/189 | 213 |
| 01_crop04 | 649/**264** | 649/264 | 670/265 | 642/263 | 670/265 | 642/263 | 264 |
| 01_crop05 | 309/**128** | 301/127 | 663/263 | 658/259 | 435/194 | 457/200 | 128 |
| 01_crop06 | 568/**265** | 554/263 | 518/240 | 504/231 | 655/263 | 661/265 | 265 |
| 01_crop07 | 632/**260** | 644/259 | 651/265 | 664/261 | 628/265 | 625/264 | 260 |

`zyx` matches 12/12 exactly. The next best, `zxy`, matches 2/12 (crop04 only, a symmetric location). This agrees with the original `batos_make_paired_crop.py` (`arr[z0:z1, y0:y1, x0:x1]`). In `crops_sibgrapi.json`, `crop_convention.axis_order = "zyx"` and `slicing_status = "verified"`.

### Crops cut (`scripts/make_crops.py`)
- 12 × `C:\BA_TOS_DATA\crops\<crop_id>\gray.tif` (raw uint16, 200³) and `labels.tif` (uint16, 200³).
- Manifest `data/metadata/crops_manifest.csv`: pixel and file SHA-256, shape, dtype, gray min/max, n_gt_labels, n_gt_touching_face, n_gt_small_lt1000, n_gt_valid.
- **Acceptance test PASSED:** n_gt_valid per crop = 265, 216, 272, 129, 271, 263 / 262, 213, 264, 128, 265, 260, identical to the historical counts. **Total 2808** = paper.
- Note: the valid-GT rule (≥ 1000 voxels, no face contact) matches both the legacy `batos_instance_metrics.py` and the e12 strict rule at margin 0.

### Repository
- `.gitignore`: `data/` → `data/*` + `!data/metadata/`. `data/raw` and `data/crops` stay ignored, and `*.tif` stays ignored globally (checked with `git check-ignore`).
- Committed on new branch `journal-extension` (see the git log for the hash) and pushed to origin. `main` is untouched. The camera-ready PDF in `paper/` was not committed (not on the requested list).

---

## 2026-11-04 — Session 4: Phase 1c, batos-paper-v2 implementation + development-crop run

### Housekeeping
- `paper/batos_sibgrapi_2026_camera_ready.pdf` committed (`5d25332`).
- `docs/BA_TOS_BASELINE_REPRODUCTION.md` written (frozen E-series pipeline, recovered vs reconstructed BG1–BG12, fidelity targets, one-crop command).

### Implementation (`src/batos_v2/`, branch `journal-extension`, one commit per stage)
| Commit | Module | Content |
|---|---|---|
| `76d25ae` | `preprocess.py` (+ config, pyproject, common) | exact e13 scale estimate; openrec background subtraction (BG2–5, 7); exact λ flattening |
| `2acc263` | `tos_select.py` | exact batos_background internal_dark; batos_tos_seed.py port (vectorized prefilter, same semantics) + true-dual split (BG8); spatialization (BG9 exclusive supports via top-down propagation; BG12 medial point); exact E11 markers |
| `7aa02ff` | `gate.py` | exact e23b/c gate + routing; exact E11 EDT-Otsu matched-count markers |
| `f32673c` | `watershed.py` | exact E11 fixed-mask watershed; exact batos_markers baseline (mask union only); size filter 500/off (ids kept) |
| `e6c4ace` | `evaluate.py` | seed-ignore per §III-C (ignored iff all seed voxels in excluded references), strict (e12) secondary; bipartite matching; pooled and crop-macro summaries in separate files |
| `d8221cf` | `run_crop.py` | orchestration, run.json provenance, results/all_runs.csv |

- New things added: `pyproject.toml` (editable install into the venv so `python -I -m batos_v2...` works). In the config, `size_filter.modes: [500, "off"]` must be quoted because PyYAML reads a bare `off` as `False`.
- Tests: 18 unit tests on synthetic volumes, all passing (`python -I -m unittest discover -s tests -t .`). They include a brute-force check of the exclusive-support propagation. The synthetic end-to-end run is a plumbing test only (24 markers for 64 synthetic grains; not tuned).
- Before running: on crop04, the r = 22 opening-by-reconstruction background is nearly constant (p5 = p50 = p95 = 15312). So BG4 `reconstructed_image` would give a flat image, and `subtract_background` (≈ offset + robust u8) is the only viable reading.

### Run: EFRGP01_00_crop04 (development crop only)
Command: `C:\BA_TOS_DATA\envs\batos-py312\Scripts\python.exe -I -m batos_v2.run_crop --crop-id EFRGP01_00_crop04 --config configs/batos_paper_v2.yaml` at commit `d8221cf` (clean tree). Output in `C:\BA_TOS_DATA\runs\batos_paper_v2\EFRGP01_00_crop04\` (run.json). Runtime 96 s, peak memory 3.1 GB.

Stage log:
- R̂ = 10.20, so openrec r = 15. For comparison, the valid-GT median volume of 12,297 voxels corresponds to an equivalent radius of 14.3.
- Boundary: dark threshold 110; internal_dark has 74,132 voxels in 20,312 components.
- ToS: 9,331,331 nodes, built in 67 s. MIN-like candidates scored ≥ τs: **1,938**. MAX-like: 455,896. Selected 1,800 MIN + 1,800 MAX.
- Markers: 256 consolidated → **255 after NMS** (the paper's per-crop range is 53–405).
- Gate: m_v = 13, q_v = 27 → **not triggered** (correct; paper's closest non-triggered crop is 17/43).

**Results (seed-ignore, primary):**
| method | size filter | F1@0.25 | F1@0.50 | F1@0.75 | P@0.50 | R@0.50 | n_pred_eval / valid GT |
|---|---|---|---|---|---|---|---|
| batos_tos (= batos_sel, gate false) | 500 | 0.734 | **0.559** | **0.491** | 0.721 | 0.456 | 172 / 272 |
| batos_tos (= batos_sel) | off | 0.734 | 0.559 | 0.491 | 0.721 | 0.456 | 172 / 272 |
| edt_otsu control | 500 / off | 0.859 | 0.697 | 0.672 | 0.800 | 0.618 | 210 / 272 |

Strict protocol: batos_tos F1@0.50 = 0.496; edt_otsu = 0.657. The size filter changes nothing on this crop: no evaluated basin is below 500 voxels, and only ignored ones differ (68 vs 83 ignored).

**Target: F1@0.50 ≈ 0.985, F1@0.75 ≈ 0.978 → NOT MET (−0.43 / −0.49).** No parameter was changed after seeing these numbers.

### Where the reconstruction diverges (read-only diagnostic on the run outputs; nothing tuned)
1. **Mask, relief and watershed are not the problem.** Frozen mask vs GT foreground: precision 0.988, recall 0.947 (mask fraction 0.647; GT foreground 0.675).
2. **The markers are clean but too few.** All 255 ToS markers lie in grains (172 on valid GT, 83 on excluded border/small GT, 0 on background), none duplicated. But only **172 of 272 valid grains get a marker**, so recall is capped at 0.63. Unseeded grains are absorbed into neighbours' basins, giving merges and IoU < 0.5, which costs precision too.
3. **The cause is the ToS MIN candidate pool.** Only 1,938 MIN-like nodes pass scoring. 62 % of the 1,800 kept have their point on background (inter-grain dark), and they reach only 172 distinct valid grains. The internal_dark mask is not limiting: it touches all 272 valid grains.
4. **This contradicts the paper's K sweep.** F1 rises from K = 3600 to 7200 (0.985 → 0.998), so the original MIN pool must have exceeded 3,600 eligible nodes (≥ 1,800 per half). Ours would saturate at 1,938. Our F1 of 0.559 sits in the paper's K = 900–1800 band (0.510–0.833), the signature of candidate starvation.
5. Ranked suspects (all in reconstructed or unverifiable steps):
   - (a) **Selector image (BG3–5, BG6 interplay):** with B ≈ constant, our selector ≈ a λ = 1 high-pass of the raw crop. If the lost openrec script produced a materially different image, the ToS changes completely. This is the main suspect.
   - (b) **True-dual scoring (BG8):** the lost script may compute MIN statistics (avg area, max height, dark percentile) differently from the legacy script, or apply K before the score threshold.
   - (c) **Representative point (BG12):** 62 % of MIN points fall on background. A different point rule (e.g. a deepest-voxel point) would change internal-component eligibility.
   - (d) **Scale estimate:** R̂ = 10.2 vs ~14.3 (r = 15 instead of ~21). Low impact here, because B is near-constant at either radius.
6. Partial evidence of correct plumbing: the EDT-Otsu control under the same mask reaches 0.697 here. Paper Table II's Otsu+WS crop-macro mean is 0.830, so the shared mask/EDT path is also somewhat below the paper on this crop, consistent with the openrec/flat image (suspect a) differing too.

Status: Phase 1c implementation complete. Fidelity is **not** reached on the development crop. The other 11 crops were not run.

---

## 2026-11-04 — Session 5: controlled reconstruction experiments on EFRGP01_00_crop04 (authorised)

Rules applied:
- Development crop only.
- K, τs, NMS, boundary rule, gate thresholds, mask and watershed are frozen. `scripts/run_v2_control.py` refuses overrides of those sections.
- One factor at a time; every hypothesis was named before its run.
- Outputs: `C:\BA_TOS_DATA\runs\batos_paper_v2_controls\<name>\result.json` (command, git state, overrides, pool, markers, seeded grains, gate, metrics).

### Controls
- **C1 (port check):** the unmodified `batos_tos_seed.py` from `origin/batos-ls-auto-tile-pilot` was run (git show → scratchpad, `-I`) on the v2 selector image with the E11 flags (`--max-selected-nodes 3600`, etc.). It gives 1,938 candidates; tree: 9,331,331 nodes, 473,167 MIN-like, max_min_height 10, avg_min_area 10.475, dark_alt_max 177. Diffed against `batos_v2.tos_select.score_candidates` (k_split min_only): **identical node list and order; max |diff| = 0 for score, area, height, alt, parent_alt, delta, support, support_area, support_steps and every sub-score; the top-1800 sets are identical. No port bug.**
- **C2 (image check)**, same selector and downstream; only the selector image changes. The two new images were made with the unmodified `src/batos/batos_preprocess.py` from main (verified identical to main) on the raw crop.
  - C2a, v2 selector (harness re-run): pool 1,938 · markers 255 · seeded 172/272 · F1@.50 0.559 · F1@.75 0.491. This reproduces the main run.
  - C2b, λ15 on raw (May image): pool 858 · markers **1** · seeded 0 · F1 0.000 / 0.000. The gate fires (1/1).
  - C2c, λ1 on raw, no openrec: pool 1,945 · markers 261 · seeded 178 · F1 0.609 / 0.538.
- **C3 (8-bit stretch):** the recovered chain before the ToS (`batos_preprocess` robust_u8) uses 0.5/99.5, which v2 already matches. The E11 p1–p99 stretch applies to the mask/control image, which v2 already ports exactly. The only unrecovered 8-bit step is the openrec `_8u` output; C3 sets it to p1–p99. Result: pool 1,902 · markers 250 · seeded 167 · F1 0.542 / 0.501. No improvement.

### Hypotheses on the opening-by-reconstruction step (named before running)
| ID | Hypothesis | Pool | Markers | Seeded | F1@.50 | F1@.75 |
|---|---|---|---|---|---|---|
| H1 | top-hat: I_f = max(I − B, 0), no median add-back | 1,906 | 256 | 176 | 0.571 | 0.527 |
| H2 | 3D Gaussian σ = 1 pre-smoothing before the opening | 1,632 | 111 | 65 | 0.137 | 0.125 |
| H3 | opening on the robust-u8 crop instead of raw uint16 | 1,885 | 261 | 173 | 0.575 | 0.512 |
| H4 | 2D per-slice opening by reconstruction, disk r | 2,421 | 261 | 168 | 0.532 (500) / 0.527 (off) | 0.431 / 0.427 |
| H5 | division I/B·median(B) instead of subtraction | 1,938 | 255 | 172 | 0.559 | 0.491 |

The EDT-Otsu control under the same mask moved with the image: 0.195 (H2) to 0.727 (H3).

### Conclusion and freeze
- Hypotheses are exhausted; **none comes within 0.01 of 0.985** (best of all variants: C2c 0.609; best openrec variant: H3 0.575).
- Every openrec variant leaves the MIN pool at ~1.6k–2.4k and seeded grains at ≤ 178/272. The paper's K sweep implies an original pool well above 3,600, so **the openrec step is not the source of the gap.**
- `configs/batos_paper_v2.yaml` is frozen with the a-priori BG defaults (subtract_background, raw, 3D, no pre-smoothing, 0.5/99.5). The new fields are written out explicitly. I did not pick the dev-crop maximum among near-equal variants.
- Evidence and the remaining suspects (true-dual MIN scoring/capacity, representative-point rule, e10 eligibility/projection) are in `docs/BATOS_PAPER_V2_SPEC.md` §3c.
- Code changes this session:
  - `preprocess.py`: new reconstructed-step options `mode` tophat/divide_background, `input` robust_u8, `dimensionality` 2d_per_slice, `presmooth_sigma`. Defaults unchanged.
  - `scripts/run_v2_control.py`: the control harness.
  - A test covering all openrec variants (19 tests pass).
- No other crop was run.

---

## 2026-11-04 — Session 6: round 2 (selector and spatialisation), S1

Authorised scope: crop04 only. Frozen parameters untouched; named hypotheses; one factor at a time.

### S1: the recovered E11 builder is the frozen marker stage (adopted verbatim)
Evidence given by the author: the E13 README lists the frozen E11 settings ("K=3600, MIN candidates, score threshold=0.35, point-internal consolidation, NMS=6, EDT weight=0.20, connectivity 1"). The E11 README states that `e11_build_markers_ablation.py` in MIN mode reproduces the E10 markers voxel by voxel.

- **Verbatim vs port.** Both unmodified scripts were run from `C:\BA_TOS_DATA\linux_import\linux_patches\BA_TOS_E11_pre_article_tests\audit\` on the crop04 main-run inputs.
  - `e11_build_markers_ablation.py --polarity min --nms-radius 6 --edt-priority-weight 0.20 --min-priority-bonus 0.05 --expected-marker <port markers>`: `regression_equal_binary_mask: True`, labels identical, 256 → 255 markers.
  - `e11_build_control_markers.py --mode otsu_edt --nms-radius 6 --match-count-from <ToS markers>`: labels identical to the port (255).
- **What the builder needs.** Inputs are `--diagnostic` (CSV) and `--internal-source` (internal_dark TIFF); it never reads an image. Required CSV columns (lines 206–213): `node, polarity, score, point_z, point_y, point_x, has_reconstructed_region`. No other column is read.
- **Point-internal consolidation:**
  - L264: point → 6-connected internal_dark component label.
  - L272: keep only polarity-matching rows whose point component is > 0.
  - L309–317: sort by (component, priority ↓, score_norm ↓, row) and keep one per component.
- **EDT weight 0.20 (L287–302):** priority = minmax(score over all valid rows, both polarities) + 0.20 · log1p(max EDT of the point's internal component) / log1p(global max component EDT) + 0.05 · [polarity == min]. Then greedy NMS r = 6 (L327).
- **`projection_distance` is hard-coded to 0.0 (L299).** No projection exists in the frozen builder.
- **Pipeline change.** New module `src/batos_v2/recovered.py` calls both scripts as subprocesses (`-I`, unmodified; script SHA-256 recorded in tos_select.json / route.json). Config: `markers.implementation: e11_verbatim`, `gate.edt_implementation: e11_verbatim`, `paths.e11_audit_dir`. The port still runs and is stored as `markers_port*` with a `port_identical_to_verbatim` flag.
- **Crop04 end-to-end check** (`S1_verbatim_markers_pipeline`): pool 1,938 · markers 255 · seeded 172/272 · F1@.50 0.559 · F1@.75 0.491 · gate 13/27. Both `port_identical_to_verbatim` flags are True.
- **Consequence:** reproduction now reduces to producing the candidate CSV that `batos_tos_seed_true_dual.py` + `e10_spatialize_candidates_no_gt.py` produced: `node, polarity, score, point_z/y/x, has_reconstructed_region` (+ `region_voxels` for the gate). S2–S5 are pending author review of this S1 report.

### S0: audit of the `--internal-source` mask on crop04 (read-only diagnostic; nothing tuned)
- **Production in v2:** `batos_v2.tos_select.internal_dark_mask`, an exact port of `batos_background.py` (compute_external_dark), on the selector image `EFRGP01_00_crop04_flatten_lambda1p0_8u.tif`. That image is the λ1 flattening of the openrec r15 image. Parameters as `run_batos_pipeline.sh` step 3: `--dark-percentile 40 --otsu-offset -10 --dark-threshold-max 135 --dark-connectivity 1`. Threshold = min(Otsu 120 − 10, P40 = 144, 135) = **110**. Dark voxels 1,400,845; border-connected (external) 1,326,713; **internal 74,132**.
- **Cross-check with main:** the unmodified `src/batos/batos_background.py` (identical to `main`) was run on the same selector, with the openrec image as `original` and the C1 seed (the seed does not affect the mask), into `runs\batos_paper_v2_controls\S0_main_batos_background`. Output: the same threshold 110 and the **voxel-identical** internal mask.
- **E13 expectations:** the mask is non-empty (the E13 validation rule passes). The selector name `*_flatten_lambda1p0_8u.tif` matches the λ1 folder convention, and no `lambda15p0` input is involved.
- **Structure (6-connected, as the E11 builder labels it):**
  - **20,312 components.** Median size 1 voxel; 14,979 are single voxels; p90 is 3 voxels; the largest is 142.
  - Dominant GT label: valid grain 8,756 · excluded grain 4,884 · background 6,672.
  - Components touching more than one valid grain: **13** (more than one GT label of any kind: 29).
  - **Every one of the 272 valid grains contains components lying entirely inside it** (8,731 such components). Median per grain 31, min 15, max 77. Grains with no component voxel: **0**.
- **Conclusion: the internal mask is not the limit.** Under one-marker-per-component, the structural ceiling is all 272 grains, with ≥ 15 candidate components each. Merging is negligible (13 shared components).
- **Where the ceiling actually comes from:**
  - Of the 1,800 selected MIN candidates, only **668 have their point inside an internal component**. They hit only **256 distinct components**, 173 of them dominated by valid grains.
  - 8,583 grain-interior components receive no candidate point.
  - The loss is in the ToS candidate points (S2–S4), not in the internal mask.
- **Side observation:** the internal mask is mostly single-voxel speckle from the λ1 high-pass. A MIN point must coincide exactly with one of these specks to survive point-internal consolidation. That makes S4 (projection of the point into the nearest internal component) directly relevant.

### H6 / H7 / closing (preprocessing, authorised before S2–S5)
- **H6** (`H6_openrec_filter_replaces_image`, openrec.mode=reconstructed_image, recovered radius r = max(8, round(1.5·R̂)) = 15; R̂ = 10.2 from the median of 1,500 EDT peaks; the p95 path is only a fallback for < 10 peaks):
  - The opening by reconstruction is ~constant (p50 = p95 = 18689) and erases every grain.
  - Internal mask: 97,622 voxels, 13,106 components, median size 1.
  - Pool 3,136 · markers 215 · seeded **0** · F1 0.000 / 0.000. **Rejected.**
- **H7:** `min_size` (`min_predicted_size`) is consumed only by `run_batos_pipeline.sh` arg 4 → `MIN_LABEL_SIZE` → `filter_predicted_labels_by_size.py --min-size` (legacy watershed output). It is not used in the background/internal mask or the candidates. `minimum_peak = max(2, P75)` is internal to the scale estimator (already exact). For crop04, min_size = max(500, round(0.1·4/3·π·10.2³)) = 500, which equals the existing 500 filter: F1 0.559 / 0.491. **H6+H7** = H6 with the 500 filter: 0.000.
- **Closing-by-reconstruction:** no hint in the recovered text or flags (only binary mask closing elsewhere). Not run.

### S2–S5 (crop04 still below target, so continued)
- **S2** (`S2_k_before_tau`, tos.k_order=k_first): pool 1,800 · markers 255 · seeded 172 · F1 0.559 / 0.491. Identical, because ≥ 1,800 MIN pass τs = 0.35. The MIN score distribution (263,774 scored) was ≥ 0.25: 77,901 · ≥ 0.30: 34,899 · ≥ 0.35: 1,938 · ≥ 0.40: 33 · ≥ 0.45: 13. So τs binds under either order, unlike the paper's flat Table IV.
- **S3** (two independent true-dual pools): already the v2 implementation (darker/brighter than parent, separate stats, each capped at 1,800). Nothing new to run. Height check: not saturated. **MIN heights are negative** (median −12, max 10), h_score = 0 for 99.99 %, support found for ~0 % of MIN nodes. The top-1800 weighted contributions are h 0.003 · delta 0.036 · area 0.154 · dark 0.168 · support 0.000.
- **Root cause (S3b):**
  - In higra 0.6.13 the 3D ToS root altitude is 0 with `padding='mean'`, same as `'zero'` (identical trees, 3,865 root children). With `'none'` it is 174. The selector border mean is 145.75.
  - `attribute_height(..., 'auto')` tests `alt_root > min(non-leaf alt)` (0 > 0 → False), so it uses "decreasing".
  - `tos.height_direction=mean_padding_emulated` reproduces the auto test with a true mean-padding root (145.75 > 0 → increasing).
  - Run `S3b_height_direction_mean_padding`: pool **204,769** · max_min_height 140 · markers 256 · seeded 172 · F1 0.563 / 0.496 · gate 12.5 / 27.
  - **With S3b the top-1800 MIN set is identical for every τs in 0.25–0.45** (≥ 0.45: 93,607), which reproduces the paper's τs plateau. Adopted in the frozen config.
- **S4** (`S4_projection_nearest_internal`, tos.point_projection=nearest_internal; points not in internal_dark are moved to the nearest internal voxel; `projection_distance` recorded): pool 1,938 · **markers 847** · seeded **223** · F1 0.658 (500) / 0.333 (off) · F1@.75 0.603 · gate 12 / 31.
  - **S3b+S4:** pool 204,769 · markers 859 · seeded 222 · F1 0.660 / 0.330.
  - Not adopted: the marker count is outside the paper's 53–405, the result depends on the size filter, and the frozen E11 builder hard-codes projection_distance = 0.
- **S5** (columns): the E11 builder reads `node, polarity, score, point_z/y/x, has_reconstructed_region`; the gate reads `region_voxels`. Every other column we write is descriptive and read by nothing in the recovered code: sub-scores, `region_touches_volume_face`, `region_internal_*`, `region_full_voxels`, and now `projection_distance` / `original_point_*`. No placeholder affects markers or routing.

### Frozen config after round 2
- `tos.height_direction: mean_padding_emulated` (adopted), `tos.k_order: threshold_first`, `tos.point_projection: none`.
- `markers.implementation` and `gate.edt_implementation`: `e11_verbatim`.
- **crop04 = 0.563 / 0.496 vs the 0.985 / 0.978 target: not reached.**
- Evidence and remaining suspects: `docs/BATOS_PAPER_V2_SPEC.md` §3d. The selector image (λ1 high-pass → 1-voxel-speckle internal mask, inconsistent with the compact markers of paper Fig. 2) and/or a different Higra version's tree are the open questions. No other crop was run.

---

## 2026-11-04 — Session 7: round 4, frozen E-series code run verbatim — REPRODUCED

Sources: `C:\BA_TOS_DATA\linux_import\chat_recovered\` (README read) and `zsh_recovered\` (INDEX read).

### Steps
1. **Frozen pipeline** (`src/batos_v2/frozen.py`, commit `6a6222f`).
   - Unmodified scripts via subprocess `-I`, flags and order from `e13_stage2_worker_20260725_165544.sh` → `run_e13_e11_all_crops.sh`. The recovered runner differs from the linux_patches copy only by `find -L` for the GT.
   - Preprocessing: e13_estimate_scale → openrec v2 (zsh_recovered) with r = max(8, round(1.5·R̂)) → run_batos_pipeline steps 1–4 with λ 1 (flags identical to the recovered Linux `run_batos_pipeline.sh`). Steps 5–6 (legacy watershed) are not consumed by E11.
   - Then true-dual K3600 → e10 spatialize → e11 markers (min, dual) and EDT-Otsu control → e11 fixed-mask watershed → e09b → e23c route + selective watershed + e09b.
   - Not run: `audit_true_batos_e07.py` (GT capacity audit, not recovered), internal_medial / mintree_e14 controls.
   - Script SHA-256 is recorded per step. The pilot-branch `batos_tos_seed.py` (needed for the background step's `--seed`) was staged unmodified to `C:\BA_TOS_DATA\frozen_code\pilot_branch_413920f\`.
2. **Golden check:**
   - The golden CSV is not crop04 (point overlap 0.03 %) and not either original crop05 (0.03 %).
   - Best match is the holdout EFRGP01_01_crop04: 29.8 % identical points, no row identity. Most likely that crop under E10-era preprocessing (unproven).
   - **Our crop05 runs vs the recovered original E13 outputs (both scans): identical** selector CSV (21 columns, 3,600 rows, same order), spatial diagnostics (41 columns), markers_tos_min and markers_edt_otsu tables, fixed-mask and surface SHA-256, and all e09b metrics (max |ΔF1| = 0).
3. **crop04 (development):** seed_ignore tos_e11_min **F1@0.50 = 0.985185, F1@0.75 = 0.977778** (targets 0.985 / 0.978). Gate off: m_v = 43, q_v = 49, with 398 markers; the E23B README states crop04 median 43. MIN ≡ DUAL. Strict protocol: 0.828. EDT-Otsu: 0.887. Runtime 369 s.
4. **Since crop04 matched:**
   - EFRGP01_00_crop05: gate **fires (2, 3)**, 53 markers. ToS F1 0.012121 / 0.000. Selective (EDT) **0.259887** / 0.237. These are the paper's ToS min 0.012 and sel min 0.260.
   - EFRGP01_01_crop05: gate **fires (2, 3)**, 191 markers. ToS 0.065574 / 0.049. Selective (EDT) 0.807273 / 0.793.
   - Holdout EFRGP01_01_crop04: gate off (43, 49), 397 markers, ToS **0.984848 / 0.984848**.
5. **Environment:** the venv is aligned to `requirements_environment_initial.txt` (numpy 2.5.1, scipy 1.18.0, scikit-image 0.26.0, higra 0.6.13, pandas 3.0.5, tifffile 2026.7.14, …). All pins match (commit `bd3ae9b`). Remaining difference: Python 3.12.10 here vs 3.13 on Linux (cpython-313 caches). With the crop05 outputs bit-identical, this has no observed effect.

### Interpretation
- The fidelity gap of rounds 1–3 came from three reconstruction errors, now replaced by recovered code:
  1. the openrec input/output stretch (8-bit p1–p99 before the opening);
  2. the exact ball;
  3. the true-dual selector itself (polarity-aware |height|, intensity score, hierarchy suppression; the legacy selector we ported is sign-sensitive to higra's direction flip).
- BG10 is answered: the E-series numbers are **without** a prediction size filter.
- The remaining 8 crops (00/01 crop02, 03, 06, 07) were not run, per instruction.

---

## 2026-11-04 — Session 8: Phase 1 closure — Table II reproduced on all 12 crops

### Runs
- **Primary runs.** The remaining 8 crops (00/01 crop02, 03, 06, 07) were run with `python -I -m batos_v2.frozen` into `C:\BA_TOS_DATA\runs\batos_paper_v2_frozen`, all rc 0, 317–455 s per crop. Code commit `6aafcc8`, clean.
- **Determinism.** 00/01 crop04 and 00/01 crop05 were rerun into `..._frozen_rerun`. e09b CSVs and ToS/sel segmentations are byte-identical to the originals.
- **Memory bug.** `frozen.py`'s first peak-memory capture measured only the venv `python.exe` launcher (about 4.7 MB per stage; on Windows the interpreter is a child process). Fixed in `ab1982f` (max `peak_wset` over the process tree, recorded as `peak_memory_method`).
- **Memory reruns.** All 12 crops were rerun into `..._frozen_mem` for memory only. Results are byte-identical to the primary runs on all 12, so the pipeline is deterministic.
  - Peak memory is 3.31–4.78 GB per crop, driven by openrec_v2; ToS/spatialize take about 3.1 GB.
  - `git.dirty = True` in those records is due only to uncommitted doc edits (CLAUDE.md, the reproduction doc) and the aggregator's `--relabel` fix. No pipeline code was modified.
- **Runtime source.** `results/all_runs.csv` takes runtime from the primary runs (sequential, uncontended) and memory from the memory runs.

### Aggregation (`scripts/aggregate_sibgrapi.py` → `results/batos_sibgrapi/`)
- **Method.** Frozen e09b unmodified. IoU 0.25 is added by `scripts/e09b_with_iou025.py`, which only overrides the default threshold tuple; its 0.50/0.75 output equals the plain e09b run on all 12 crops (`checks.csv`). Primary protocol: seed_ignore.
- **Table II (crop-macro, no size filter, sample SD).** All statistics are within ±0.0005:
  - ToS 0.8238 ± 0.3670, median 0.9850, min 0.0121, F1@0.75 0.8175 (paper 0.824 ± 0.367 / 0.985 / 0.012 / 0.817).
  - sel 0.9063 ± 0.2099, median 0.9850, min 0.2599, F1@0.75 0.8992 (paper 0.906 ± 0.210 / 0.985 / 0.260 / 0.899).
  - Otsu+WS 0.8301 ± 0.1819, median 0.8816, min 0.2599, F1@0.75 0.7698 (paper 0.830 ± 0.182 / 0.882 / 0.260 / 0.770).
  - The paper's SD is the sample SD; population SDs are 0.008–0.016 lower.
- **Gate.** It fires on 00_crop05 and 01_crop05 only, both (2, 3). The closest non-triggered crop is 01_crop03 at (17, 43). ToS and sel segmentations are bitwise identical on the 10 gate-negative crops.
- **ToS markers per crop.** 53–405.
- **Pooled object (seed_ignore, F1@0.50).**

  | Route | No filter | Filter 500 (labelled) |
  |---|---|---|
  | sel | **0.9496** (TP 2630, FP 101, FN 178, P 0.963, R 0.937) | 0.9522 (FP 86) |
  | ToS | 0.9116 (P 0.932, R 0.892) | 0.9351 (FP 45) |
  | Otsu+WS | 0.8657 | 0.8657 (no object < 500) |

  - The legacy 0.974 (May pipeline) is not comparable; it is reported alongside, not as a target.
- **Over-segmentation ratio (pred_eval / valid GT, pooled, no filter).** sel 0.973, ToS 0.957, Otsu+WS 1.095.
- **Size filter.** Without `--relabel`, the filter removes only FPs (TP unchanged everywhere).
  - Table II matches **filter off** (filtered ToS mean 0.8297, min 0.0153), so the published E-series numbers did not use it (BG10 confirmed on 12 crops).
  - A first aggregation with `--relabel` was **discarded**: e09b links predictions to seeds by label value, and relabelling corrupted the seed-ignore sets (ToS crop04 0.726 with 4 objects removed).
- **Not compared.** Min-tree-core+WS (0.768, E14 control not run), Table IV sweeps, Table V crop extent.

### Bookkeeping
- `results/all_runs.csv` now holds the 72 frozen rows: 12 crops × 3 routes × filter on/off, with runtime and peak memory. The previous Phase 1c port rows were moved to `results/all_runs_phase1c_port.csv`.
- CLAUDE.md §2 (filter row, Phase 1 status) and `docs/BA_TOS_BASELINE_REPRODUCTION.md` were updated.
- Tag `batos-paper-v2-frozen`.

---

## 2026-11-05 — Session 9: Phase 1d, label-agnostic evaluator (`batos_v2.evaluate_v2`)

### Why
Frozen e09b decides seed-ignore by label value: it labels the marker volume (6-conn), classifies each marker
component by its dominant GT label, then drops the *prediction label with the same value*. This breaks for any
method without markers (Cellpose, SAM, nnInteractive) and for any relabelling (observed in Session 8: `--relabel`
took ToS crop04 from 0.985 to 0.726).

### What (`src/batos_v2/evaluate_v2.py`, new; no frozen code touched)
- **Reference protocol.** Identical to e09b: valid GT ≥ 1000 voxels and not touching any face; excluded = the rest; IoU computed on the domain `~ignore_mask`.
- **Support region S(p).**
  - With markers: the marker voxels inside p (geometric seed region; label values are irrelevant).
  - Without markers: p's own voxels.
- **Ignore rules.**
  - `dominant` (default with markers): the most frequent GT label among the GT-positive voxels of S(p) is excluded. Ties go to the smallest label, as in e09b.
  - `fraction` (default without markers): |S(p) ∩ excluded| / |S(p)| ≥ `--ignore-overlap-threshold` (default 0.50).
  - Background-supported predictions are never ignored, so they stay FP. Unseeded predictions (no marker voxel) are evaluated.
- **Matching.** Maximum cardinality at IoU ≥ t (TP identical to e09b); among maximum matchings, the one with maximum total IoU (`linear_sum_assignment` per connected component). This makes mean matched IoU label-invariant. e09b keeps scipy's arbitrary maximum matching, so its `mean_matched_iou025` depends on label order. At t ≥ 0.5 the matching is unique.
- **Added metrics.**
  - IoU 0.25.
  - `ap_iouXXX` = TP/(TP+FP+FN): Cellpose/StarDist convention, no scores needed.
  - `coco_ap_iouXXX`: 101-point, score-ranked; NaN unless `--scores` is given. The plan does not define AP; both are recorded and labelled.
  - Mean matched IoU.
  - `pred_to_gt_ratio`.
  - Morphology per object matched at IoU 0.50, on full objects and in voxel units unless `--voxel-size` is given: volume, d_eq, centroid, surface area (marching cubes), sphericity, elongation (sqrt of the largest/smallest covariance eigenvalue), each with its error (pred − GT).
- **Tests.** `tests/test_evaluate_v2.py` (12 synthetic tests, incl. equivalence to the frozen e09b functions, label invariance, background-seeded FP, rule thresholds, matching, COCO AP and morphology on a ball/ellipsoid). Full suite 31/31 OK.

### Verification (`scripts/verify_evaluate_v2.py` → `results/evaluate_v2/`)
- **Reference.** The Phase 1 closure e09b outputs (frozen e09b plus the IoU 0.25 wrapper) on 12 crops × {BA-TOS_ToS, BA-TOS_sel, Otsu+WS} × size filter {none, 500} = 72 evaluations.

| Check | Pass |
|---|---|
| A: markers + dominant ≡ e09b seed_ignore (counts, P/R/F1 @ .25/.50/.75, mean matched IoU @ .50/.75) | **72/72** |
| B: prediction labels randomly permuted (+ offset) ≡ A, all columns incl. mIoU@.25 | **72/72** |
| C: no markers, fraction 0.50 ≡ e09b overlap_ignore | **72/72** |
| D: markers + fraction 0.50 ≡ e09b seed_ignore (the two rules coincide on these data) | **72/72** |

- **mIoU@0.25 vs e09b.** v2 − e09b ranges 0 to +0.0112 (nonzero in 62/72 runs). This is expected from the max-IoU tie-break; excluded from A by design.
- **Before the matching fix.** B failed 29/72, only on mIoU@0.25. This confirms that e09b's mIoU@0.25 is label-order dependent.

### What the evaluator gives without markers (own-voxel fraction 0.50; = e09b overlap_ignore)

| Route (no filter) | Crop-macro F1@0.50, seed (markers) | Crop-macro F1@0.50, no markers | Pooled F1@0.50, seed / no markers |
|---|---|---|---|
| BA-TOS_ToS | 0.8238 ± 0.3670 | 0.8162 ± 0.3616 | 0.9116 / 0.9021 |
| BA-TOS_sel | 0.9063 ± 0.2099 | 0.8947 ± 0.2197 (min 0.2130) | 0.9496 / 0.9401 |
| Otsu+WS | 0.8301 ± 0.1819 | 0.8120 ± 0.1902 | 0.8657 / 0.8484 |

- **Effect of dropping markers.** Without markers, scores are about 0.01–0.02 lower. A prediction seeded inside a valid grain but lying ≥ 50 % over excluded GT gets ignored, and its valid grain becomes a FN (e.g., pooled sel: TP 2630 → 2574, FP 101 → 94).
- **Comparing methods.** Methods without markers (Cellpose etc.) must be compared against BA-TOS evaluated in the **same** no-marker mode, not against the seed-ignore numbers.

### Morphology (matched @0.50, pooled, voxel units, no filter, BA-TOS_sel)
| Measure | Median |
|---|---|
| Volume relative error | −0.052 (abs 0.060) |
| d_eq error | −0.51 vox (GT median d_eq 28.7 vox) |
| Centroid error | 0.08 vox (p95 8.5) |
| Sphericity error | −0.027 (GT 0.870) |
| Elongation error | +0.006 |

- Predictions are systematically ~5 % smaller than GT. ToS and Otsu+WS are similar (Otsu+WS centroid p95 10.4).

### Bookkeeping
- Per-object morphology (`morphology_matched_seed.csv.gz`, ~7 MB) is gitignored. The tracked summaries are `morphology_summary_pooled.csv`, `per_crop_all_modes.csv`, `crop_macro_by_mode.csv`, `pooled_object_by_mode.csv`, `verification.csv` and `REPORT.md`.
- No rows were added to `results/all_runs.csv`: this is a re-evaluation of existing runs, not a new method run.

## 2026-11-05 — Phase 5a group 1 (selector): power outage and recovery

- Group 1 launched 13:26:02 (`python -I -m batos_v2.sensitivity --identity --groups selector --workers 10 --ram-gb 22`,
  408 crop-jobs = 12 identity + 33 perturbations × 12; `_git.json`: 8cceddf + uncommitted Phase 5a code, committed
  unchanged as 4a4f259 at 13:28; the YAML was edited at 13:26:32, which does not affect identity jobs). A power outage
  stopped it after 8 jobs (last at 13:43:04).
- Completeness criterion: `metrics.csv` (evaluate_v2) + `run.json` with end `timestamp` (+ `DONE`). Kept: identity on
  00_crop03/04/05/07, 01_crop02/03/04/05. Deleted as partial: identity 00_crop02 (stopped at spatialize),
  00_crop06 (watershed), 01_crop06 (flatten), 01_crop07 (openrec); `selector_K_1800` 00_crop02..07 (≤ true_dual).
  No perturbed job had finished. The 4 `gate_median_limit_*` jobs (group-4 smoke test, 13:25) are complete and kept.
- Identity resumed on the 4 missing crops at 13:55 (commit 4a4f259, same budget, ~9.5 min). **Identity verified 12/12:**
  per-crop F1@0.50/0.75 identical to `results/evaluate_v2/per_crop_all_modes.csv` (max |Δ| = 0); sel 0.9063 ± 0.2099
  (min 0.2599), ToS 0.8238 ± 0.3670 (median 0.985, min 0.0121), F1@0.75 0.8175 → 0.8992; gate only on both crop05
  (2, 3), closest off 01_crop03 (17, 43) == Table II.
- Group 1 perturbations (396 crop-jobs) relaunched with the same settings; log `C:\BA_TOS_DATA\runs\sensitivity\group1_selector.log`.

### 2026-11-05 16:45 — group 1 killed by Claude Code memory-pressure reaper; relaunched at 16 GB / 8 workers
- Killed at 16:37 (system near 32 GB while the session was idle; not a pipeline failure). 212/396 perturbation
  crop-jobs complete, 0 failed. Deleted 10 partials: `selector_area_max_4400` {00_crop04, 01_crop04, 01_crop06,
  01_crop07}, `selector_dark_percentile_max_39` 00_crop02..07. Relaunched `--groups selector --workers 8 --ram-gb 16`
  (user decision; deterministic, so only throughput changes). Group 2 will use the same 16 GB / 8 workers.

### Gate blind spot found in group 1 (for Phase 5b)
Under `selector_w_area` 0.18→0.36 and `selector_w_delta` 0.12→0.06 the gate stops firing on **EFRGP01_00_crop05**
(01_crop05 still fires, mv=2, qv=3):

| run | 00_crop05 gate | mv | qv | ToS markers / preds | sel F1@0.50 | ToS F1@0.50 |
|---|---|---|---|---|---|---|
| identity | fires | 2 | 3 | 53 preds | 0.2599 | 0.0121 |
| w_area 0.36 | silent | 18 | 18 | 1 pred | 0.0000 | 0.0000 |
| w_delta 0.06 | silent | 14.5 | 16.25 | 2 preds (1 ignored) | 0.0000 | 0.0000 |

(129 valid GT.) The ToS stage collapses to 1–2 markers; median/Q75 of their region voxels are then taken over 1–2
values and come out above the 3/4 limits, so the gate does not fire. The exclusive-support gate tests support *size*
but not marker *count*, so a count collapse gets past it. Closest non-triggered crop elsewhere stays far away (01_crop03:
44/56 and 39/54). A Phase 5b gate should also test the marker count (e.g. against the GT-free R̂ expected count).

### 2026-11-05 19:33 — group 1 (selector) complete
- 408/408 crop-jobs complete, 0 failed (resumed run: 396 jobs at 16 GB / 8 workers, system peak 12.3 GB above start,
  absolute 26.2 GB). Ranked with `scripts/rank_sensitivity.py` → `results/sensitivity/{global_sensitivity.csv,
  ranking.csv, RANKING.md}`; report `results/sensitivity/GROUP1_SELECTOR_REPORT.md`.
- 14/17 selector parameters influential (> 0.01). Top: area_preferred_position 0.613, w_delta 0.512, w_area 0.417,
  area_preferred_width 0.417, area_min 0.289, height_scale_percentile 0.216, K 0.127. Not influential: area_factor,
  min_score, w_support (0.000, dead-code check confirmed).
- Gate: beyond the two cases above, `delta_norm 40` (mv 11, qv 14.5, 3 preds) and `area_min 3` (mv 4, qv 7, 21 preds)
  also silence it on 00_crop05 (sel F1 0); `area_preferred_position 0.84` makes it fire on 01_crop03 (mv 2, qv 4;
  sel 0.652 vs ToS 0.542). The `area_min 3` case shows that the 3/4 limits are tied to `area_min = 2`.
- Sensitivity runs are not appended to `results/all_runs.csv` (provenance is per-job `run.json` + the collected CSV);
  open question for the user.
- Next: group 2 (openrec_flatten) at 16 GB / 8 workers.

### 2026-11-05 22:31 — group 2 (openrec_flatten): main run done; radius_factor x2 infeasible on this machine
- Main run 19:34–22:28 (16 GB / 8 workers): 121/132 ok. 10/11 perturbations complete on 12 crops. System peak
  18.7 GB above start, absolute 31.9 GB (= all RAM; openrec is the heavy stage).
- `openrec_radius_factor = 3.0` (radius 24–31 instead of 12–16): 11 crops failed with `MemoryError` in the frozen
  openrec script (`opening_by_reconstruction`, line 35, scipy grey erosion with a ball footprint). Only 01_crop03
  (r = 24) finished (openrec 732 s).
- Sequential re-run (`--workers 1`, 22:28): r = 27–31 fail with MemoryError within ~1 s each, alone, so the failure is
  deterministic and does not depend on concurrency. 01_crop05 (r = 24) was running when Claude Code's memory-pressure
  guard killed the shell; its partial folder was deleted. Status: 1/12 done, 10 FAILED (r 27–31), 01_crop05 not run.
- Groups 3 and 4 not started; waiting for the user's decision on how to treat radius_factor x2.

### 2026-11-05 23:58 — radius_factor 2.0 extra point killed by the memory-pressure guard (stopped, per user: stop on failure)
- `--values 2.0 --workers 2 --ram-gb 16` (commit cd947a8), started 23:02. 9/12 complete (00_crop02..07, 01_crop02,
  01_crop03, 01_crop05; r 16–20). System peak 17.8 GB above start with only 2 workers. Killed at ~23:58 while
  01_crop04 (r 20, at step 7/12, past openrec) and 01_crop06 (r 21, in openrec) were running; 01_crop07 not started.
  Partial folders 01_crop04 and 01_crop06 left in place (the driver deletes them on relaunch). Group 2 not ranked or
  committed; groups 3–4 not started.

### 2026-11-06 — group 2 (openrec_flatten) complete and ranked; overnight terminal run interrupted at 09:02
- The user ran the remaining sweeps in a separate terminal (2026-11-06 00:30). radius_factor 2.0: 12/12 (the 3 missing
  crops with `--workers 1`, done 00:53). Group 2: 133/144 complete; radius_factor 3.0 is infeasible (10 MemoryError at
  r ≥ 27, 01_crop05 at r = 24 not run).
- Ranking: 5/6 influential. radius_factor 0.883 (0.75 → all crops at r = 8, mean 0.024), lambda_xy 0.637 (λ 0.5),
  radius_floor 0.036, out_percentile 0.026, sigma_z 0.016; norm_percentile 0.004 not influential. Gate: fires on 10/12 at
  radius_factor 0.75, but the fallback cannot rescue a preprocessing collapse; silent on 00_crop05 at λ 2.0 (4/7) and
  σz 1 (3/5). Report: `results/sensitivity/GROUP2_OPENREC_FLATTEN_REPORT.md`.
- Group 3 (markers) in the terminal: 682/1500 crop-jobs done by 09:02 (all `legacy_*`, bg_dark_percentile 20, and 10/12
  of 80), then interrupted (exit 0xC000013A = Ctrl+C / console closed). bg_dark_percentile 80 on both crop05 is a real
  break: e11 watershed "marker volume empty after applying the mask" (rc 1). bg_otsu_offset −20 / −5: 13 jobs
  FAILED by the interrupt and 8 partial; the driver re-runs these on relaunch. Group 4 not started.

## 2026-11-06 — Phase 5c: extinction-based ToS selector on EFRGP01_00_crop04 (development crop only)
- Code 790b0a6: `batos_v2.extinction_selector` replaces only the frozen true_dual stage; `batos_v2.extinction_run`
  reuses sensitivity.CropRun (frozen upstream outputs, frozen downstream calls, evaluate_v2 with markers / seed-ignore).
  Command: `python -I -m batos_v2.extinction_run --crop EFRGP01_00_crop04` (variants volume, dynamics, product; K 3600).
  Ran alongside the user's group 3–4 terminal sweep; this work peaked at 3.13 GB per stage (driver 0.62 GB).
- F1@0.50 (ToS = sel, gate off for all): volume 0.906, product 0.764, dynamics 0.611; frozen 0.985. Target ≥ 0.98 not met.
  Cause: split grains (47 / 109 / 192 valid grains with ≥ 2 markers vs 0 frozen); recall 0.96–0.99.
  Pool 238,168 dark extrema (0 touch a face); e10 accepted every node (alt / parent-alt checks).
- Parameters: 1 numeric (K) + variant vs 17 in the frozen selector.
- Side check (GT for diagnosis only): e13 r_est 10.20 / 9.00 / 8.00 (crop04 / 00_crop05 / 01_crop05) vs median GT
  r_eq 14.32 / 14.35 / 14.36 → underestimates by 29 / 37 / 44 %, more on the gate-firing crops.
- Write-up: `docs/EXTINCTION_SELECTOR.md`; tables `results/extinction/`; 6 rows appended to `results/all_runs.csv`.

## 2026-11-06 — Phase 5c follow-up A (uniqueness / scale window) and Phase 5b side study B (scale estimator)
- Code 37c9b83 (+ driver failure handling, this commit). Crops: 00_crop04, 00_crop05, 01_crop05 only. Every stage
  peaked ≤ 3.66 GB, run in sequence. Logs: `C:\BA_TOS_DATA\runs\followup_AB.log`, `followup_A2.log`.
- A (crop04 F1@0.50): volume_unique 0.460 (1,307 markers, 244 split: 19,230 internal fragments), scale window
  r_est 0.065 / r_ps_peak 0.111 (39–57 markers), literal grey-volume window 0.508. Frozen 0.985, volume v1 0.906.
  Target 0.98 not met. volume_unique_scale breaks on both crop05 (0 candidates → e10 KeyError on the empty CSV;
  recorded FAILED, deterministic).
- B: pattern-spectrum peak of e13's Otsu foreground = 10 / 11 / 11; e13 r_est 10.20 / 9.00 / 8.00. GT median
  **inscribed** radius 10.63 on all three crops (GT r_eq 14.3; non-spherical grains, ratio 1.35). This corrects the
  earlier side check: e13 underestimates on crop05 only (−15 / −25 %). EDT local-max median 2.0–2.2 (noise).
- B rerun (openrec radius 15 / 16 / 16 via `openrec_lowmem`, bit-identical to frozen at r = 15): crop04 identical to frozen
  (0.985). crop05 gate still fires (2 / 4, 2 / 3), ToS 0.007 / 0.016. sel 0.722 (from 0.260) / 0.702 (from 0.807),
  driven by the count-matched EDT-Otsu fallback (ToS markers 53 → 254, 191 → 364). The scale estimate is NOT the root
  cause of the crop05 collapse.
- Docs: `docs/SCALE_ESTIMATOR.md` (new), `docs/EXTINCTION_SELECTOR.md` §7. Tables: `results/scale/`,
  `results/extinction/summary.csv`; 18 rows upserted in `results/all_runs.csv`.

## 2026-11-06 — B2 crop05 forensics (read-only)
- The user's message was truncated after "(b) … Name the first"; this was read as "the first stage where crop05 diverges".
  A second announced run was not included in the message → not done.
- (a) The 16 % unlabelled bright foreground is the specimen's confining membrane: one 26-connected sheet (99.9 %),
  ≤ ~6 voxels thick, at raw grey ≈ 31.5k (grains 47.3k, air ≈ 9k). It runs diagonally and cuts off the high-y /
  high-x corner of both crop05 (same origin, consecutive scans). About 45 % of each crop05 is exterior air. Renders:
  `C:\BA_TOS_DATA\runs\crop05_forensics\*_midslices.png`.
- (b) First divergence = flatten λ = 1 (`batos_preprocess.flatten_volume_u8`): bg_ref = median(blur[blur > 0]) = 0.9 on
  crop05 vs 186 on crop04. The openrec image is 58–67 % zero, and the float Gaussian tails make 51 % of the "positive"
  voxels < 1. The selector loses its offset (grain body 44 / 33 vs 165; cores / pores / air 0–7). Downstream: dark_alt_max
  46 / 42 vs 178; MIN points 97 % / 86 % in GT background; internal mask = grain-body fragments (83 / 75 per grain);
  53 / 191 markers with 2-voxel supports → gate. Verified by re-running the frozen flatten (bit-identical).
  Doc: `docs/CROP05_FORENSICS.md`; table: `results/crop05_forensics/stage_comparison.csv`.

## 2026-11-06 — Root-cause fix test (flatten bg_ref) and A2 (extinction + spatial NMS)
- Code 1f4b84b. Crops: 00_crop04, 00_crop05, 01_crop05. Stages peaked ≤ 3.36 GB, run in sequence. Log:
  `C:\BA_TOS_DATA\runs\rootcause_fix.log`. Frozen code untouched.
- Fix: copy of the λ=1 flatten (`batos_v2.flatten_fixed`), bg_ref over the specimen. A = median(blur ≥ 1), B = median
  over e13's Otsu foreground of the raw crop. Re-run from flatten onward, all later stages frozen.
  - bg_ref A / B: crop04 198.8 / 223.8; 00_crop05 176.7 / 210.6; 01_crop05 117.4 / 162.0 (frozen 186.1 / 0.9 / 0.9).
  - crop04 unchanged: 398 markers, gate 43/49, F1@0.50 0.985, F1@0.75 0.978.
  - 00_crop05: 193 markers, gate off (43/49), ToS = sel 0.992 / 0.992 (frozen ToS 0.012, sel 0.260).
  - 01_crop05: 196 markers, gate off (8/21), ToS = sel 0.913 / 0.891 (frozen ToS 0.066, sel 0.807).
  - A and B give identical results: the robust rescale cancels the offset shift.
  - **Conclusion: the crop05 support collapse is a flattening artefact caused by exterior air, not a property of the method.**
  - Hybrid indicative macro (crop05 replaced only): ToS 0.976 ± 0.023. Not a v2 result; the 12-crop rerun is pending.
- A2 volume_nms on crop04 (r = 7 / 10 / 13):
  - markers 494 / 536 / 499; seeded 266 / 263 / 260; split 61 / 53 / 10.
  - F1@0.50 0.860 / 0.803 / 0.804, against volume without NMS 0.906 and frozen 0.985. NMS does not help.
- Doc: `docs/CROP05_ROOT_CAUSE.md`. Tables: `results/flatten_fix/summary.csv`, `results/extinction/summary.csv`, `results/all_runs.csv`.

## 2026-11-06 — batos-paper-v2-fixflat on all 12 crops
- Frozen batos-paper-v2 with one change: the λ=1 selector flatten offset is bg_ref = median(blur ≥ 1)
  (`batos_v2.flatten_fixed --variant A`). Re-run from the flatten onward; scale and openrec are read from the frozen
  run, every later stage is frozen.
  - Driver: `flatten_fix_run` (gained `--root/--label/--summary/--no-all-runs`).
  - Run sequentially, max stage peak 3.36 GB, alongside the user's Phase 5a sweep (groups 3–4 not touched).
  - Runs: `C:\BA_TOS_DATA\runs\batos_paper_v2_fixflat\A`. Log: `C:\BA_TOS_DATA\runs\fixflat_12.log`.
- The gate fires on none of the 12 crops; ToS = sel bitwise everywhere.
  - Both crop05: 43/49 and 8/21, against frozen 2/3.
  - Other 10: gate stats and F1 at all three IoUs identical to frozen; ToS segmentation bitwise identical on 9 of 10.
  - 01_crop03 changes 0.9 % of voxels with the same TP/FP/FN.
- Crop-macro F1@0.50, filter none: ToS = sel 0.976 ± 0.023, median 0.986, min 0.913 (01_crop05); F1@0.75 0.970.
  - Frozen: ToS 0.824 ± 0.367, sel 0.906 ± 0.210.
  - Filter 500: 0.983 ± 0.014, min 0.957.
  - crop05: 00 0.992 (frozen ToS 0.012 / sel 0.260); 01 0.913 (0.066 / 0.807).
- Pooled F1@0.50 (2808 GT):
  - filter none: 0.979 (P 0.979, R 0.978, TP 2746, FP 58, FN 62); frozen ToS 0.912 / sel 0.950.
  - filter 500: 0.984 (P 0.990, FP 29, FN 62); frozen ToS 0.935 / sel 0.952. Legacy 0.974 is a different pipeline.
- Runtime: 87.6 min for 12 crops, against frozen 71.7.
  - Scan 00 crops ran under the sweep's CPU contention; scan 01 crops took 242–268 s from the flatten onward,
    close to frozen. The flatten itself takes about 2 s.
- Otsu+WS control not re-run (not produced by this driver; unchanged).
- Development-set result; to be frozen before any test crop.
- Report: `results/batos_v2_fixflat/REPORT.md`. Scripts: `scripts/aggregate_fixflat.py`, `scripts/check_fixflat_identity.py`.
  48 rows in `results/all_runs.csv`. Tag `batos-paper-v2-fixflat`.

## 2026-11-06/07 — Phase 1f: test protocol, frozen vs fixflat on 52 held-out crops (6 specimens)

- Method under test frozen beforehand: `configs/batos_paper_v2_fixflat.yaml` (tag `batos-paper-v2-fixflat`). No parameter
  changes. Inventory `data/metadata/dataset_inventory.md`; crops `data/metadata/crops_test.json` (GT-free rule, seed
  20261106, committed 4a72daa before cutting); cut by `scripts/extract_and_cut.py` (52 crops, 8,614 valid GT; extracts
  deleted, EFR `gt_full` kept).
- Run: `python -I scripts/run_test_protocol.py` (sequential, peak ≤ 3.7 GB), 16:03 → 23:42. All 51 runs used the frozen
  openrec_v2 (lowmem never triggered). e13 padded to 12 rows per specimen, padding checked identical.
- Deviation: the first pipeline crash (fixflat EMSGP01_00_test01, e11 "no MIN candidate inside internal_dark") stopped the
  run at crop 29; the driver now records crashes in `FAILED.json` and scores them as empty segmentations (F1 0), then
  resumed. 5 crashes total (4 fixflat, 1 frozen: EOSGP01_00_test03, where fixflat could not run either; scored 0).
- Results (crop-macro F1@0.50, sel, filter off): EFR_internal frozen 0.879 ± 0.223 → fixflat 0.956 ± 0.048 (filter 500:
  0.985 ± 0.021); EFS 0.882 → 0.968; EMR 0.303 → 0.250; EMS 0.432 → 0.479; EOR 0.293 → 0.238; EOS 0.265 → 0.365;
  all 52: 0.538 → 0.575. Pooled F1@0.50: EFR 0.929 → 0.967; all 0.700 → 0.752. Fixflat gate never fires on EF; frozen
  gate fires on 5 EF boundary crops (crop05 signature reproduced).
- Forensics (`scripts/forensics_test_protocol.py`, `results/test_protocol/forensics_stages.csv`): 25 flagged crops. The
  M/O collapse diverges at stage 03 (frozen `batos_background.py`): dark threshold min(Otsu−10, p40, 135) falls to the
  p40 branch, the dark phase percolates, all seeds rejected as external. Branch separates fixflat F1 completely
  (p40: 20 crops ≤ 0.21; Otsu: 31 crops ≥ 0.65). For interior crops the branch follows the openrec radius (15 → Otsu,
  ≤ 12 → p40), and e13 underestimates the radius ≈ 2× on M/O grains (7–8 vs GT inscribed 14.5–15.3; granulometric
  estimate 14–15). The fixflat gain is the same stage: 14 of 21 frozen-p40 boundary crops move to the Otsu branch.
- Nothing tuned. Report `results/test_protocol/REPORT.md`; 832 rows in `results/all_runs.csv`
  (experiment_id `test_protocol_<method>/<crop>`). `results/sensitivity/global_sensitivity.csv` left uncommitted (user sweep).
- Next (needs the user's decision): Phase 5d (queued); any scale/stage-03 remedy is a method change → dev design, re-freeze,
  new held-out set.

## 2026-11-07 — Phase 5e step 1: v3 scale source r_granulo (dev crops only)

- Decision (user): BA-TOS v3 (scale-free) designed on the 12 EFR dev crops + a synthetic scale test only; the 52 test
  crops stay untouched until v3 is frozen.
- `src/batos_v3/scale.py`: two-pass r_granulo = argmax of the binary pattern spectrum (exact EDT ball openings,
  r = 2..40) of the specimen-only Otsu foreground; specimen = fill_holes(LCC(closing(fg0, ball(2 r0)))).
- Iteration 0 (raw Otsu foreground): r_granulo = 7 on all 12 dev crops (enclosed noise cavities cap the EDT).
  Iteration 1: fill_holes (parameter-free) on both passes → r_granulo = 10 on all 12 dev crops (r0 = 10–11), openrec
  radius 15 (= v2's), specimen fraction crop05 0.49 / crop03 0.78 / crop06 0.97 (cf. geometry estimate 0.50/0.79/0.92).
  r_ref (dev median) = 10.
- openrec_lowmem timing at r = 22 on a dev crop: 157 s (kept for r > 15).

## 2026-11-07 — Phase 5e step 2: v3 constants, patched copies, driver (dev crops only)

- `src/batos_v3/constants.py` (every pipeline constant, kind L/V/N/VF/D/Q/R/X/T, v3 rule; R_REF = 10),
  `src/batos_v3/make_patched.py` → `src/batos_v3/patched/` (e11 watershed + controls: min_component 500 and closing
  iterations 1 exposed; e23c: gate limits 3/4 exposed; defaults = frozen), `src/batos_v3/pipeline.py` (presets
  fixflat and v3). Design record: `docs/V3_SCALE_FREE.md`.
- Intensity-guard audit over 115 dev+test runs: stage-03 center-dark-max / min-ring-contrast 0 rejections (removed);
  stage-04 center-dark-max 3 rejections of 1.93 M (removed); stage-04 min-ring-contrast 130k rejections (kept, R).
- Design iteration 2 (before any v3 run): gate limits scaled v2·max(1, s³) instead of v2·s³ (at s = 0.7 v2·s³ gives
  m_v ≤ 1.03, which can no longer detect a 2–3-voxel collapse).
- Checks: fixflat preset bit-identical to batos-paper-v2-fixflat on EFRGP01_00_crop04 (ToS and sel SHA-256).
  v3 on EFRGP01_00_crop04: F1@0.50 0.985 (stage-03 threshold 117 vs v2 107); EFRGP01_00_crop05: 1.000 (115 vs 97).
- Synthetic crops: `scripts/make_synthetic_scale.py` → 36 crops (12 dev × f 0.7 / 1.4 / 2.0, centred affine resample
  of the full volumes, linear gray / nearest labels; f = 1 identity verified against every dev crop),
  `data/metadata/crops_synth_scale.json`, `crops_manifest_synth.csv`; valid GT ≥ round(1000 f³).

## 2026-11-07 — Phase 5e step 3: synthetic scale test PAUSED (checkpoint, partial)

- Paused at 14:2x on user request (PC shutdown). Run killed cleanly between crops (EFRGP01_00_crop03_x2 was in openrec; its
  partial run dir will be wiped and redone on resume). 25 of 48 crop evaluations complete: f = 1.0 (12), f = 1.4 (12), f = 2.0 (1).
- Partial crop-macro F1@0.50 (markers/dominant, no filter, BA-TOS_sel; ToS identical):
  f 1.0 fixflat 0.976 ± 0.023 / v3 0.989 ± 0.008; f 1.4 fixflat 0.750 ± 0.372 (2 crashes "no MIN candidate", scored 0) / v3 0.995 ± 0.007;
  f 2.0 (1 crop, 00_crop02) fixflat 0.000 / v3 1.000. Step-3 criterion at 1.0 met (no design iteration needed so far).
- Peak memory ≤ ~3.1 GB per stage. Tables in `results/v3_scale_test/` are PARTIAL and will be regenerated at the end.
- Resume: `C:\BA_TOS_DATA\envs\batos-py312\Scripts\python.exe -I scripts/run_v3_scale_test.py` (skips DONE/EMPTY/FAILED runs and
  existing eval JSONs), then `python -I scripts/aggregate_v3_scale_test.py`. Remaining: 11 crops at f 2.0, 12 at f 0.7.

## 2026-11-08 — Phase 5e step 3: synthetic scale test COMPLETE (48/48)

- Resumed 2026-11-07 23:53. The one partial run dir (`fixflat/x2/EFRGP01_00_crop03_x2`, stopped in openrec) was deleted;
  the 25 complete evaluations were kept; the 23 missing crops (11 × f 2.0, 12 × f 0.7) were rerun with the same settings,
  END 05:01. Ran concurrently with the user's sensitivity sweep (`--groups markers`, 6 workers), so runtimes are inflated.
- Crop-macro F1@0.50 sel, no filter, fixflat → v3: f 0.7 0.999 → 0.975; f 1.0 0.976 → 0.989; f 1.4 0.750 → 0.995;
  f 2.0 0.239 → 0.938. Own size filter (500 / 500 s³): 1.000 → 0.995; 0.983 → 0.990; 0.824 → 1.000; 0.383 → 1.000.
  Pooled F1@0.50 (no filter): 0.999 → 0.978; 0.979 → 0.988; 0.819 → 0.995; 0.277 → 0.947.
- r_granulo = 7 / 10 / 14 / 20–21 vs GT r_in 7.35 / 10.63 / 14.44 / 20.52; e13 saturates (10.0 at f 1.4, 12.6 at f 2.0).
  fixflat crashes ("no MIN candidate inside internal_dark"): 2 at f 1.4, 6 at f 2.0. v3: 0 crashes, 0 empty. Gate never fires.
- v3 below 0.90 (no filter): both crop05 at f 2.0 (0.667, 0.621; 9 valid GT, 9 and 11 sub-filter fragments; with the
  500 s³ filter 1.000). v3 is behind fixflat at f 0.7 without the filter (351 small FP fragments vs 13; filtered 0.995 vs 1.000).
  Reported, not tuned. Step-3 criterion met; no design iteration. Report `results/v3_scale_test/REPORT.md`.

## 2026-11-08 — Phase 5e step 4: frozen v3 on the 52 held-out test crops (COMPLETE)

- `python -I scripts/run_test_v3.py`; freeze hash check passed in every session (code commit a853e83, tag
  batos-paper-v3-scalefree). Four sessions: machine shutdown after 10 crops (partial EFRGP01_01_test11 deleted);
  stopped at 14:36 by Claude Code's low-memory stop after 27 (partial EMRGP01_01_test08); user-terminal session from 17:19, stopped after 29
  (partial EMSGP01_00_test02); final session 20:43 → END 23:14. Partials rerun from scratch, completed crops never re-run.
  `test_v3.log` mixes UTF-8/UTF-16 lines (PowerShell redirection in session 3).
- 52/52 DONE, 0 crashes, 0 empty, gate never fires (frozen 30, fixflat 13). Aggregated with the frozen
  `scripts/aggregate_test_v3.py`, unmodified; 624 rows `test_v3_v3/<crop>` in `results/all_runs.csv`.
- Crop-macro F1@0.50 sel (seed-ignore, no filter) frozen / fixflat / v3: all 0.538 / 0.575 / **0.912**; EFR_internal
  0.879 / 0.956 / 0.985; external 0.435 / 0.460 / 0.890; EMR 0.980, EMS 0.986, EFS 0.974, EOR 0.754, EOS 0.754 (v3).
  Pooled F1@0.50: 0.700 / 0.752 / **0.921** (filter 500: 0.722 / 0.760 / 0.922).
- r_granulo / GT r_in 0.94–0.98 on all six specimens (Pearson 0.998, n = 52, three size clusters); e13 0.46–0.95 (r −0.81).
- v3 runtime median 397 s/crop (frozen 271), 6.0 h total, peak 3.65 GB; openrec_lowmem on all 32 EM/EO crops. Not contention-controlled.
- 10 flags < 0.90: 8 EO interior + EOSGP01_01_test08 (boundary; only v3 < fixflat case, 0.776 vs 0.982) + EFSGP01_01_test05
  (0.876, identical for all three methods, crop-intrinsic). The frozen first-stage rule says "watershed" for all. Supplementary
  `scripts/forensics_test_v3_seed_density.py` (diagnosis only) puts the first anomaly at **candidate capacity**: MIN
  candidates are at the 1800/s³ cap on all 52 crops; EO interior has 1.19 candidates / 0.35 markers per grain (valid +
  excluded GT) vs 2.7 / 0.60 on EFR. Grain counts per crop fall as ~0.70 (not s⁻³ = 0.30) because face-touching grains
  (~60 %) scale as a surface. Not acted on: a K-scaling change is a method change → dev design, re-freeze, new held-out set.
- Report `results/test_v3/REPORT.md`. `results/sensitivity/global_sensitivity.csv` still left uncommitted (user sweep).

## 2026-11-09 — Phase 5f steps 1–2: v3.1 (count constants floored at v2) — validation FAILED, stopped before freeze

- Step 1: `src/batos_v31` changes only kind N in v3: v2 · max(1, s⁻³). This is the user's choice over c · N̂, which
  scales as s⁻³ like v3. N constants: tos_seed max-selected-nodes 900, markers_legacy max-markers 420, true_dual
  max-selected-nodes 3600 (K = 1800 MIN + 1800 MAX). Implemented by rebinding `batos_v3.pipeline.constants`. No
  frozen file was edited.
- Step 2 (`scripts/run_v31_dev.py`, `scripts/aggregate_v31_dev.py`, results/v31_dev/):
  - Flags are identical to v3 for all 24 crops with s ≤ 1, and crop04 re-ran bit-identical. f = 1.0 stays 0.989 (pass).
  - Crop-macro sel F1@0.50 (no filter): f = 1.4 0.995 → 0.941; f = 2.0 0.938 → 0.326 (12/12 crops < 0.90).
  - crop05 at f = 2.0: 0.667 / 0.621 → 0.247 / 0.290. **Fail.**
  - Cause: markers per grain 0.42 → 2.62 at f = 2.0, giving over-segmentation into sub-filter fragments. FN = 0
    everywhere, and with the own filter F1 is 0.996.
- Steps 3–4 not run: no freeze, no tag, test set 2 not defined or cut (still unseen). Decision needed from the user.

## 2026-11-09 — Decision: v3 is the journal method; Phase 5a sweep complete and ranked

- Decision (user): `batos-paper-v3-scalefree` is the journal's method. The EO capacity limit goes in the limitations
  section with the seed-density analysis (`results/test_v3/REPORT.md`). v3.2 (per-grain merge/suppression) is deferred.
  v3.1 stays unfrozen; test set 2 stays unseen.
- Phase 5a: the user's sweep finished. Groups 3 (markers, 1497/1500 crop-jobs) and 4 (gate_watershed, 216/228) are
  complete; every configured perturbation was attempted. `scripts/rank_sensitivity.py` re-run on the complete
  CSV: 96/96 method parameters ranked, 25 influential. Summary: `results/sensitivity/FULL_RANKING.md`.
- Failures: `bg_dark_percentile` 80 (both crop05) and `bg_dark_connectivity` 3 (01_crop05) are real breaks
  (e11 empty markers). `ws_mask_closing_iterations` 0 failed on all 12 crops. This is a sweep artifact: scipy
  `iterations=0` means "until convergence" and empties the mask. The point is invalid and was not re-run.
- Next: Phase 2 (Cellpose baselines).

## 2026-11-09/10 — Phase 5h: v3.2 per-grain marker uniqueness (worktree `..\batos-v32`, branch `journal-extension-v32`)

- Design record and acceptance criteria were written before any run: `docs/V32_UNIQUENESS.md`.
  - Hypotheses are tested in the order H1 → H2 → H3; the first that holds is frozen.
  - Thresholds: f = 1.0 crop-macro sel F1@0.50 ≥ 0.985; f = 0.7 / 1.4 / 2.0 ≥ v3 − 0.005.
- Code: `src/batos_v32/` (uniqueness rules, supports export, prefix/tail split), `scripts/run_v32_dev.py`,
  `scripts/aggregate_v32_dev.py`.
  - The shared v3.1 prefix runs once per crop and is kept unpruned for later hypotheses.
  - H1 basin seeds = the EDT-Otsu control's markers before count matching (NMS 6 s). No new constant.
- Identity: with hypothesis `none` the prefix/tail split reproduces v3 (EFRGP01_00_crop04, f = 1) and v3.1
  (EFRGP01_00_crop02_x2, f = 2) bit-for-bit, on both segmentations.
- H1, f = 1.0 complete (12/12). sel, seed-ignore, no filter:
  - crop-macro F1@0.50 0.991 ± 0.008 (v3 0.989), min 0.974 (00_crop06, = v3); F1@0.75 0.988 (v3 0.986).
  - pooled F1@0.50 0.989: TP 2767, FP 18, FN 41 (v3: 2765 / 25 / 43).
  - Grains: valid grains split 5 → 0; unseeded 27 → 27.
  - **Criterion 1 passes.** Run continues with f = 2.0, 1.4, 0.7 (`C:\BA_TOS_DATA\runs\v32_dev_H1.log`).
- H1, f = 2.0 complete (12/12; sel, seed-ignore, no filter): crop-macro F1@0.50 **0.876** ± 0.086 (min 0.692,
  00_crop05; 7 crops < 0.90). v3 0.938, v3.1 0.326. Threshold 0.933 → **criterion 2 FAILS; H1 does not hold.**
  - Pooled: TP 188, FP 50, FN 0 (v3 188/21/0; v3.1 188/819/0).
  - All 188 valid grains are seeded. 30 are split (2.2 % → 15 % median per crop); v3 14, v3.1 184.
  - 19 markers lie in background.
  - Every FP is a sub-filter fragment: with the own filter (500 s³) H1 = v3 = 1.000.
  - Cause: the EDT-Otsu basins (seeds NMS 6 s = 0.6 r_granulo) put 2–3 basins in one large grain, so a grain keeps
    one marker per basin. Reported, not tuned (pre-registered order: next is H2).
  - H1 continues at f = 1.4 / 0.7 for the full report; those prefixes are shared with H2/H3.
- H1, f = 1.4 complete (12/12): crop-macro F1@0.50 **0.985** ± 0.027 (min 0.901, 00_crop05). v3 0.995, v3.1 0.941.
  Threshold 0.990 → criterion 2 also fails at f = 1.4.
  - Pooled: TP 877, FP 21, FN 0 (v3 877/9/0; v3.1 877/90/0). Split valid grains 6 (v3 3, v3.1 60).
  - 16 markers in background (v3 7). With the own filter: 1.000 (all FP are sub-filter fragments).
- H1, f = 0.7 complete (12/12; supplementary): crop-macro F1@0.50 0.977 (v3 0.975; threshold 0.970) → passes.
  F1@0.75 0.977 (v3 0.974). Split median 0.5 % → 0.2 %.
- **H1 verdict: does not hold** (passes f = 0.7 and 1.0; fails f = 1.4 and 2.0). Tables: `results/v32_dev/tables.md`.
  Next per pre-registered order: H2 (NMS 0.8 r_granulo on the final marker list), reusing the same prefixes.
- H2 (NMS 0.8 r_granulo on the final list), all factors 12/12. sel, seed-ignore, no filter:
  - f = 0.7: 0.977 (pass). f = 1.0: 0.991 (pass; pooled TP 2767 / FP 18 / FN 41, same as H1).
  - f = 1.4: **0.951** (min 0.755; pooled FP 72 vs v3 9) → fails. f = 2.0: **0.441** (pooled TP 188 / FP 492 / FN 0;
    3.98 markers per valid grain, 100 % of valid grains split) → fails.
  - The radius is too small for these grains: markers 0.8 r_granulo apart still fit several per grain, so H2
    removes far fewer v3.1 duplicates than H1.
  - **H2 does not hold.** Next per the pre-registered order: H3 (merge touching exclusive supports ≤ V_g).
- H3 (merge touching exclusive supports ≤ V_g), all factors 12/12. Peak 3.6 GB (supports export).
  - f = 0.7 0.975 (pass), f = 1.0 0.989 (pass), f = 1.4 **0.941** (fail), f = 2.0 **0.326** (fail).
    At every factor the result equals v3.1 bit-for-bit in score.
  - Cause: 0 merges on all 48 crops. `touching_pairs` = 0 everywhere.
  - Checked on EFRGP01_00_crop02_x2 to rule out a bug:
    - all 213 marker labels are present in the leaf volume;
    - leaf voxel counts equal `region_voxels`;
    - each marker point lies inside its own support.
  - The supports are compact cores of ~100–300 voxels, against V_g = 11 494 at f = 1.4. They cover 5.5 % of the
    crop. Only 6 of the 3 600 candidate supports touch any other, and none of those are kept markers.
  - So with exclusive supports, the rule as specified cannot act. **H3 does not hold.**
- **Phase 5h outcome: no hypothesis holds** (all three pass f = 0.7 / 1.0 and fail f = 1.4 / 2.0).
  - Per the pre-registered rule: no freeze, no `batos-paper-v3.2` tag.
  - Test set 2 was not cut and is still unseen; seed 20261109 is reserved.
  - v3 (`batos-paper-v3-scalefree`) stays the journal method.
  - Tables: `results/v32_dev/tables.md`.
