# Crop05 root cause — flatten fix test and A2 (extinction + spatial NMS)

Date 2026-11-06. Crops: EFRGP01_00_crop04, EFRGP01_00_crop05, EFRGP01_01_crop05 only. Frozen code untouched; GT used
for evaluation only. Every stage peaked ≤ 3.36 GB, run in sequence (< 6 GB budget). Log:
`C:\BA_TOS_DATA\runs\rootcause_fix.log`. Code (commit 1f4b84b): `src/batos_v2/flatten_fixed.py`,
`src/batos_v2/flatten_fix_run.py`, `src/batos_v2/extinction_selector.py` (variant `volume_nms`). Tables:
`results/flatten_fix/summary.csv`, `results/extinction/summary.csv`; rows in `results/all_runs.csv`
(`phase5_flatten_fix_{A,B}/*`, `extinction_volume_nms_r{7,10,13}`). Runs: `C:\BA_TOS_DATA\runs\flatten_fix\{A,B}\<crop>`.

## 1. Fix test: bg_ref computed over the specimen only

`flatten_fixed.py` is a copy of the frozen λ = 1 flatten (`src/batos/batos_preprocess.py`, `flatten_volume_u8`). It
differs only in the offset that is added back:

- **frozen:** `bg_ref = median(blur[blur > 0])`. The copy with `--variant frozen` reproduces the frozen selector bit for bit on all three crops.
- **A:** `bg_ref = median(blur[blur ≥ 1])`, which ignores the sub-grey-level Gaussian tails over the zero plateaux.
- **B:** `bg_ref = median(blur[Otsu foreground of the raw crop])`, using e13's own foreground (imported read-only).

Upstream stages (scale, openrec) are read from the frozen run. Every downstream stage is the unmodified frozen call:
legacy ToS seed, background, legacy markers, true_dual, e10, e11 markers, EDT-Otsu, watershed, E23c gate, selective
watershed and evaluate_v2.

| crop | variant | bg_ref (frozen) | ToS markers | grains seeded / split | gate mv / qv | gate | ToS F1@0.50 / @0.75 | sel F1@0.50 / @0.75 |
|---|---|---|---|---|---|---|---|---|
| 00_crop04 | frozen | 186.1 | 398 | 268 / 0 | 43 / 49 | off | 0.985 / 0.978 | 0.985 / 0.978 |
| 00_crop04 | **A** | 198.8 | 398 | 268 / 0 | 43 / 49 | off | **0.985 / 0.978** | 0.985 / 0.978 |
| 00_crop04 | **B** | 223.8 | 398 | 268 / 0 | 43 / 49 | off | **0.985 / 0.978** | 0.985 / 0.978 |
| 00_crop05 | frozen | 0.9 | 53 | 32 / – | 2 / 3 | **fires** | 0.012 / – | 0.260 / – |
| 00_crop05 | **A** | 176.7 | 193 | 129 / 0 | 43 / 49 | off | **0.992 / 0.992** | 0.992 / 0.992 |
| 00_crop05 | **B** | 210.6 | 193 | 129 / 0 | 43 / 49 | off | **0.992 / 0.992** | 0.992 / 0.992 |
| 01_crop05 | frozen | 0.9 | 191 | 85 / – | 2 / 3 | **fires** | 0.066 / – | 0.807 / – |
| 01_crop05 | **A** | 117.4 | 196 | 125 / 11 | 8 / 21 | off | **0.913 / 0.891** | 0.913 / 0.891 |
| 01_crop05 | **B** | 162.0 | 196 | 125 / 11 | 8 / 21 | off | **0.913 / 0.891** | 0.913 / 0.891 |

Notes:
- Frozen crop05 "grains seeded" counts the distinct valid grains hit by MIN points (`docs/CROP05_FORENSICS.md`, stage 5).
- The 43 / 49 on 00_crop05 is the real gate statistic. It equals crop04's by coincidence: the q25 differs (24 vs 19).

**crop05 recovers on the ToS route.**
- On both crop05 the gate no longer fires. The ToS markers have real supports: median 43 and 8 voxels, against 2 in the frozen run.
- ToS F1@0.50 goes from 0.012 to **0.992** on 00_crop05 and from 0.066 to **0.913** on 01_crop05. Both are above the frozen sel route (0.260 / 0.807).
- crop04 stays exactly at 0.985 / 0.978 under both variants: same 398 markers, same gate statistics.

**The support collapse on crop05 is a flattening artefact caused by exterior air, not a property of the method.**
- About 45 % of crop05 is air beyond the confining membrane, and openrec clips it to exactly 0.
- The `> 0` guard then lets the float Gaussian tails over those zero plateaux pull the median offset down to 0.9.
- Computing the offset over the specimen restores the selector, and the frozen BA-TOS chain does the rest.

**A and B are equivalent.** The offset only shifts `vol − blur`, and the robust p0.5–p99.5 rescale that follows
cancels the shift whenever the clip at 0 is inactive. The selectors from A and B are therefore identical on 00_crop05
and crop04. On 01_crop05 they differ by at most one grey level (rounding), and the downstream results are identical.
Any offset taken over the specimen works. A (`blur ≥ 1`) needs no extra input and no extra parameter.

### Consequences (indicative only, not a v2 result)
- With the fix, the E23c gate fires on **none** of the three crops. On crop05, the gate and the EDT-Otsu fallback were
  compensating for this artefact.
- If the two crop05 ToS values in Table II are replaced by the fixed ones, the crop-macro ToS F1@0.50 goes from
  0.824 ± 0.367 (min 0.012) to **0.976 ± 0.023 (min 0.913, 01_crop05)**, and F1@0.75 from 0.817 to 0.970.
- This is a hybrid figure: the other ten crops were **not** re-run with the fix. Their offsets are not near-zero
  (crop04: 186 → 199, with an unchanged result), but that has to be checked crop by crop before it is reported. Next
  step: run variant A on all 12 crops, then on the test crops after freezing.

## 2. A2: extinction (volume) + greedy spatial NMS, crop04 only

`volume_nms` is the volume-extinction selector from `docs/EXTINCTION_SELECTOR.md`, with two additions:
- **Spatial NMS.** Before the top-K cut, the ranked dark extrema are suppressed greedily: an extremum is dropped if
  its flat-zone centroid lies within radius r of an already kept one. Selection stops at K = 3600.
- **Representative choice.** The face rejection and the representative rule are unchanged.

Everything downstream is frozen. Radius 10 is the granulometric r (`r_ps_peak`, crop04); 7 and 13 bracket it.

| selector | NMS r | NMS survivors (≤ K) | ToS markers | markers in valid / bg / excluded | grains seeded (of 272) | split grains | gate mv / qv | P / R @0.50 | F1@0.50 | F1@0.75 |
|---|---|---|---|---|---|---|---|---|---|---|
| frozen true_dual | – | – | 398 | – | 268 | 0 | 43 / 49 | – | **0.985** | 0.978 |
| volume (no NMS) | – | 3600 | 477 | 323 / 1 / 153 | 271 | 47 | 45 / 168 | 0.833 / 0.993 | 0.906 | 0.906 |
| volume_nms | 7 | 3600 | 494 | 335 / 2 / 157 | 266 | 61 | 157.5 / 184 | 0.777 / 0.963 | **0.860** | 0.851 |
| volume_nms | 10 | 3600 | 536 | 320 / 36 / 180 | 263 | 53 | 153 / 182 | 0.708 / 0.926 | 0.803 | 0.793 |
| volume_nms | 13 | 2420 | 499 | 270 / 85 / 144 | 260 | 10 | 158 / 186 | 0.710 / 0.926 | 0.804 | 0.794 |

**NMS does not help. Every radius is below the plain volume selector (0.906) and well below frozen (0.985).**
- Spatial NMS spreads the selection more evenly. That pushes MIN points into the pores and background (36–85 background
  markers at r = 10–13, against 1 without NMS). It costs precision, and it costs recall too (260–266 grains seeded, against 271 without NMS).
- Split grains drop only at r = 13 (10), but the radius is then larger than the inscribed radius (10.6). Neighbouring
  grains suppress each other: 2420 survivors < K, and 12 grains are unseeded.
- The splits do not come from duplicate extrema close together. They come from several internal-dark components per
  grain (one marker each, e11), which a point-distance NMS on extremum centroids cannot merge.
- A2 therefore does not close the gap to 0.985. Volume extinction remains a weaker selector than the frozen weighted
  score on crop04.

## 3. Summary
1. The crop05 failure is an implementation artefact of the selector flatten (bg_ref over the "positive" blur, polluted by
   exterior-air zero plateaux). A one-line, GT-free fix (offset over `blur ≥ 1` or over the specimen foreground)
   recovers ToS F1@0.50 0.992 / 0.913 on the two crop05 without changing crop04 (0.985).
2. With the fix, the support-collapse gate is not triggered on these crops. Its role in the paper has to be re-assessed
   once the fix has been run on all 12 crops.
3. Extinction + spatial NMS (r = 7 / 10 / 13) gives 0.860 / 0.803 / 0.804 on crop04 and is not a replacement for the frozen
   selector.
