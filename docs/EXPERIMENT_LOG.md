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
