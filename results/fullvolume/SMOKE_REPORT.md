# Phase 5g step 4 — smoke test `smoke_center3` (3×3×3 central tiles), frozen v3

Volume EFRGP01_00 (1600×1400×1400 zyx). Tiles t{2..4}_{2..4}_{2..4}: cores of 200³ covering z,y,x 400–1000. Each core
gets a halo of 20 (2·r_granulo, with r_granulo = 10 estimated GT-free on the central 200³). That gives a padded tile
of 240³. The method is frozen `batos-paper-v3-scalefree` (hash-checked), BA-TOS_sel route, and the gate fired in
0/27 tiles. Valid GT: the label is in `labelsInside.txt`, has ≥ 1000 voxels, and its bbox lies inside the 600³
core box. That leaves **6,005** GT (of 9,278 ROI labels in the whole volume). Evaluation is the streaming
evaluate_v2 (seed-ignore, dominant rule), verified identical to evaluate_v2 on 3 dev crops.

## Pooled object metrics (no size filter; filter 500 in metrics.csv)

| stitch variant | F1@.25 | F1@.50 | F1@.75 | P@.50 | R@.50 | pred eval | GT in merges | GT split | seam F1 | interior F1 |
|---|---|---|---|---|---|---|---|---|---|---|
| spec (centroid ownership + IoU ≥ 0.5 one-to-one, core authority) | 0.900 | 0.827 | 0.800 | 0.913 | 0.756 | 4,970 | 16.2 % | 0.2 % | 0.727 | 0.868 |
| + orphan rescue | 0.938 | 0.865 | 0.831 | 0.915 | 0.820 | 5,381 | 16.9 % | 0.2 % | 0.858 | 0.869 |

Seam check (spec variant): 28.8 % of matched GT were matched by a reconciled prediction, and 31 % of valid GT
cross a core boundary.

## Finding 1: seam grains owned by no tile (fixed by orphan rescue)

Among GT whose centroid is within 5 voxels of a core plane, the FN rate is 64.8 % and **49.6 % have no
overlapping prediction at all**. Further than 20 voxels from a plane the FN rate is 18 % and the no-overlap rate
is 1 %. Forensics on 4 such grains shows that both adjacent tiles segment the grain well (≈ 11.6k of ≈ 12.3k
voxels). However, each tile's instance is biased toward its own halo, so its centroid falls in the *other*
tile's core and neither tile owns it. In tile t2_2_3 (core y < 600) the centroid is at cy 601.5; in t2_3_3 it is
at cy 596.8. This bias is systematic, so centroid ownership alone is not reliable at seams.

Orphan rescue is an optional stitch step, deterministic and GT-free (`stitch --orphan-rescue`). Each non-owned
instance whose centroid lies in a selected core and whose voxels are still ≥ 50 % unlabelled is written into the
free voxels only. 719 instances were rescued. Seam F1 rises from 0.727 to 0.858, matching interior F1 (0.869).
The rescued instances use a 1-voxel proxy seed (the free voxel nearest the centroid). This has the same
cardinality as the stitched seeds (one marker voxel per instance), but it is an approximation. Tiles run from now
on store `seeds_all_<route>.npz`, so a rescued instance can use its real markers.

## Finding 2: interior recall is far below the crops (frozen-K density)

The two dev crops inside the block were evaluated with the frozen evaluate_v2 on 200³ cuts of the stitched
volume (the same crop GT and crop protocol):

| crop | v3 run on the crop | stitched tiles, same box |
|---|---|---|
| EFRGP01_00_crop04 | F1@.50 0.985 | F1 0.604, R 0.724 |
| EFRGP01_00_crop07 | F1@.50 0.987 | F1 0.646, R 0.745 |

Precision in the cut box is not comparable: stitched instances cut by the box face lose their seed and count as
FP. Recall is comparable. Even far from seams, 18 % of GT are FN, mostly with a best IoU of 0.4–0.5, which
indicates under-segmentation (16 % of GT lie in merged predictions vs 0.2 % split).

Cause: the true_dual candidate pool is a fixed count. It saturates at K = 3600 (1800 MIN + 1800 MAX) both on a
200³ crop and on a 240³ tile (`e11/spatial/spatialization_summary.json`). A 240³ tile holds 1.73× the volume, so
it gets 1/1.73 = 0.58 of the per-grain candidate density the frozen K was calibrated on. Tiles produce about
600 ToS markers (crops ≈ 395). The dev K sweep at half capacity (K = 1800) gave 0.833 on crop04; the smoke gives
0.827. This is the EO capacity limit already logged for v3 (2026-11-09). v3 scales K with s⁻³ but not with
volume.

Not acted on: changing K is a method change. A tiling that keeps the frozen regime is a **padded tile of 200³**
(core 160 + halo 20), the exact volume the method was designed and frozen on. That is a tiling change, not a
method change, and needs the user's decision.

## Resources (machine shared with 2 other heavy sessions)

27 tiles, all DONE. Wall time 13,996 s (518 s per tile with 3 workers and a stage-level memory gate). Median
time per tile 1,522 s under contention; the solo tile took 793 s. Peak RSS per tile 5.32–5.47 GB, and the
process-tree peak was 8.9 GB (cap 11 GB). Stitching took 19 s (27 s with rescue) and evaluation 8 s.
