# batos-paper-v2-fixflat — 12 crops

Date 2026-11-06. **batos-paper-v2-fixflat** is the frozen `batos-paper-v2` pipeline with **one change**. The λ = 1
selector flatten takes its offset over the specimen, `bg_ref = median(blur ≥ 1)`, instead of `median(blur > 0)`. The
change is in the copy `src/batos_v2/flatten_fixed.py --variant A`; the frozen file is untouched. Root cause and the
three-crop test: `docs/CROP05_ROOT_CAUSE.md`.

- **Inputs:** scale and openrec are read from the frozen runs.
- **Downstream:** every stage after the flatten is the unmodified frozen call. That covers legacy ToS seed, background,
  legacy markers, true_dual K = 3600, e10, e11 ToS-MIN markers, EDT-Otsu, fixed-mask watershed, E23c gate, selective
  watershed.
- **Evaluation:** the frozen e09b with the IoU 0.25 wrapper, reused from `scripts/aggregate_sibgrapi.py`. Seed-ignore
  protocol; size filter off (primary, as in Table II) and 500 (no relabel).
- **GT:** used for evaluation only. No parameter was tuned.

**How it was run**
- Runs: `C:\BA_TOS_DATA\runs\batos_paper_v2_fixflat\A\<crop>`, in sequence, alongside the user's Phase 5a sweep.
- Driver: `python -I -m batos_v2.flatten_fix_run --variants A --crops <12> --root C:\BA_TOS_DATA\runs\batos_paper_v2_fixflat --label batos_paper_v2_fixflat --summary results/batos_v2_fixflat/run_summary.csv --no-all-runs`
- Aggregation: `python -I scripts/aggregate_fixflat.py`, then `python -I scripts/check_fixflat_identity.py`.
- Log: `C:\BA_TOS_DATA\runs\fixflat_12.log`.

## Table II — crop-macro F1 (12 crops, seed-ignore, sample SD)

| route | size filter | F1@0.50 mean ± SD | median | min (crop) | F1@0.75 mean | F1@0.25 mean |
|---|---|---|---|---|---|---|
| **BA-TOS_ToS fixflat** | none | **0.976 ± 0.023** | 0.986 | 0.913 (01_crop05) | **0.970** | 0.986 |
| BA-TOS_sel fixflat | none | 0.976 ± 0.023 | 0.986 | 0.913 (01_crop05) | 0.970 | 0.986 |
| BA-TOS_ToS frozen (Table II) | none | 0.824 ± 0.367 | 0.985 | 0.012 (00_crop05) | 0.817 | – |
| BA-TOS_sel frozen (Table II) | none | 0.906 ± 0.210 | 0.985 | 0.260 (00_crop05) | 0.899 | – |
| BA-TOS_ToS fixflat | 500 | 0.983 ± 0.014 | 0.989 | 0.957 (01_crop05) | 0.977 | 0.993 |
| BA-TOS_sel fixflat | 500 | 0.983 ± 0.014 | 0.989 | 0.957 (01_crop05) | 0.977 | 0.993 |
| BA-TOS_ToS frozen | 500 | 0.830 ± 0.361 | 0.987 | 0.015 | 0.822 | – |
| BA-TOS_sel frozen | 500 | 0.909 ± 0.211 | 0.987 | 0.260 | 0.902 | – |

**The gate never fires, so ToS and sel are bitwise identical on all 12 crops.** The Otsu+WS control was not re-run:
this driver builds EDT-Otsu only as the gate fallback, so the frozen 0.830 control is unaffected and not repeated here.

## Per crop (BA-TOS_ToS = BA-TOS_sel, size filter none)

| crop | bg_ref | ToS markers | gate mv / qv (frozen) | GT | pred | TP | FP | FN | F1@0.50 | F1@0.75 | frozen ToS / sel F1@0.50 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 00_crop02 | 192.8 | 405 | 44 / 50 (44 / 50) | 265 | 266 | 263 | 3 | 2 | 0.991 | 0.987 | 0.991 / 0.991 |
| 00_crop03 | 182.8 | 299 | 42 / 49 (42 / 49) | 216 | 216 | 213 | 3 | 3 | 0.986 | 0.986 | 0.986 / 0.986 |
| 00_crop04 | 198.8 | 398 | 43 / 49 (43 / 49) | 272 | 268 | 266 | 2 | 6 | 0.985 | 0.978 | 0.985 / 0.985 |
| **00_crop05** | 176.7 | 193 | **43 / 49 (2 / 3, fired)** | 129 | 131 | 129 | 2 | 0 | **0.992** | 0.992 | 0.012 / 0.260 |
| 00_crop06 | 184.5 | 371 | 38 / 50 (38 / 50) | 271 | 267 | 260 | 7 | 11 | 0.967 | 0.967 | 0.967 / 0.967 |
| 00_crop07 | 193.2 | 393 | 43 / 50 (43 / 50) | 263 | 262 | 259 | 3 | 4 | 0.987 | 0.987 | 0.987 / 0.987 |
| 01_crop02 | 194.1 | 404 | 43 / 49 (43 / 49) | 262 | 262 | 257 | 5 | 5 | 0.981 | 0.973 | 0.981 / 0.981 |
| 01_crop03 | 170.5 | 285 | 17 / 43 (17 / 43) | 213 | 208 | 200 | 8 | 13 | 0.950 | 0.922 | 0.950 / 0.950 |
| 01_crop04 | 198.6 | 397 | 43 / 49 (43 / 49) | 264 | 264 | 260 | 4 | 4 | 0.985 | 0.985 | 0.985 / 0.985 |
| **01_crop05** | 117.4 | 196 | **8 / 21 (2 / 3, fired)** | 128 | 137 | 121 | 16 | 7 | **0.913** | 0.891 | 0.066 / 0.807 |
| 01_crop06 | 190.6 | 376 | 41 / 49 (41 / 49) | 265 | 264 | 261 | 3 | 4 | 0.987 | 0.987 | 0.987 / 0.987 |
| 01_crop07 | 194.6 | 392 | 42 / 49 (42 / 49) | 260 | 259 | 257 | 2 | 3 | 0.990 | 0.990 | 0.990 / 0.990 |


**The fix leaves the 10 non-crop05 crops unchanged.**
- **Gate statistics:** identical to frozen on all 10.
- **F1@0.25 / 0.50 / 0.75:** identical to frozen, with filter off and with filter 500.
- **Selector image:** bitwise identical on 4 of the 10 crops. On the other 6 it differs by at most 1 grey level, on
  ≤ 0.003 % of the voxels (01_crop03: 1.7 %). These are rounding effects of the shifted offset.
- **ToS segmentation:** bitwise identical to frozen on 9 of the 10. 01_crop03 changes 0.9 % of the voxels, with the
  same TP/FP/FN.
- Details: `identity_vs_frozen.csv`.

All of the gain comes from the two crop05. The closest non-triggered crop is still 01_crop03 (17 / 43). The second
closest is now 01_crop05 (8 / 21), which stays clear of the gate thresholds (≤ 3 ∧ ≤ 4).

## Pooled objects (2808 valid GT, 12 crops, seed-ignore)

| route | size filter | pred | TP | FP | FN | P@0.50 | R@0.50 | F1@0.50 | F1@0.75 | F1@0.25 |
|---|---|---|---|---|---|---|---|---|---|---|
| **ToS = sel fixflat** | none | 2804 | 2746 | 58 | 62 | 0.979 | 0.978 | **0.979** | 0.974 | 0.988 |
| **ToS = sel fixflat** | 500 | 2775 | 2746 | 29 | 62 | 0.990 | 0.978 | **0.984** | 0.979 | 0.993 |
| ToS frozen | none | 2688 | 2505 | 183 | 303 | 0.932 | 0.892 | 0.912 | 0.907 | 0.923 |
| ToS frozen | 500 | 2550 | 2505 | 45 | 303 | 0.982 | 0.892 | 0.935 | 0.930 | 0.947 |
| sel frozen | none | 2731 | 2630 | 101 | 178 | 0.963 | 0.937 | 0.950 | 0.944 | 0.969 |
| sel frozen | 500 | 2716 | 2630 | 86 | 178 | 0.968 | 0.937 | 0.952 | 0.947 | 0.972 |

The legacy May-2026 pipeline gave a pooled 0.974 (P 0.980, R 0.967, FP 39, FN 77). That is a different pipeline,
reported for reference only. fixflat with filter 500 gives 0.984 (P 0.990, R 0.978, FP 29, FN 62).

## Runtime and memory

| | fixflat (scale + openrec from frozen + flatten onward) | frozen |
|---|---|---|
| total, 12 crops | 5,254 s (87.6 min) | 4,303 s (71.7 min) |
| per crop, median | 426 s | 367 s |
| per crop, range | 280–633 s | 267–455 s |
| peak memory of the largest stage | 3.36 GB (ToS stages) | – |

- The fixflat end-to-end time is the stages from the flatten onward plus the frozen scale and openrec seconds.
- **Scan 00 crops: 428–558 s** from the flatten onward. They ran while the user's sweep (6 workers) was competing for the CPU.
- **Scan 01 crops: 242–268 s**, close to frozen.
- The fixed flatten itself takes about 2 s. **The timing difference is contention, not the fix.**
- Per-crop times: `gate_and_runtime.csv`.

## Files
- `per_crop_metrics_seed_ignore.csv` and `per_crop_metrics_all_protocols.csv`
- `crop_macro_summary.csv` and `pooled_object_summary.csv`
- `table2_fixflat_vs_frozen.csv`
- `gate_and_runtime.csv`
- `identity_vs_frozen.csv`
- `checks.csv`: e09b equals the driver's evaluate_v2 at IoU 0.50 / 0.75, 24 of 24
- `run_summary.csv`
- `results/all_runs.csv`: 48 rows, `batos_paper_v2_fixflat/<crop>` × ToS / sel × filter none / 500

## Caveat
The 12 crops are the development set that the root cause was diagnosed on. The fix has no tunable parameter and does
not use GT. Even so, these numbers are a development result: it must be frozen before any test crop is evaluated.
