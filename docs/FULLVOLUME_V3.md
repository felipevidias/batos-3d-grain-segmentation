# Full-volume tiled BA-TOS v3 — rules (Phase 5g)

The method is frozen `batos-paper-v3-scalefree`. Its files are hash-checked before every run, and none is edited.
The rules below are tiling / stitching / evaluation rules of the full-volume application. All were fixed a priori
and GT-free. Config: `configs/fullvolume_v3.yaml`. Driver: `scripts/run_fullvolume.py`. Launcher:
`scripts/run_full_c160.cmd`.

1. **Scale is a specimen property, estimated once per scan (decided 2026-11-10).**
   - r_granulo comes from `batos_v3.scale.estimate`, run GT-free on the central 200³ crop of the volume (EFRGP01_00: r = 10).
   - Every tile uses it: the frozen scale stage still runs on each padded tile, because it also writes the tile's
     specimen mask used by v3's stage-03 threshold.
   - Then `scale.json` is rewritten with r_granulo = 10 and openrec_radius = max(8, round(1.5 r)) = 15. The tile's
     own estimate is kept in `scale_tile_estimate.json` and in tile.json (`tile_r_granulo_own`).
   - Why: on tiles with little specimen the per-tile estimate is meaningless. t0_1_2 (specimen fraction 0.05) gave
     r 39 and openrec radius 58, which stalled the first launch (archived as `full_c160_aborted_pertile_scale`).
     t0_1_3 gave r 2. On the 27 specimen-filled smoke tiles it was always 10.
   - Note: the per-tile specimen mask is still built by the frozen stage with the tile's own r0.
2. **Halo** = ceil(2 · r_granulo) = 20.
3. **Tile geometry (full run):** cores of 160³ plus the halo, i.e. 200³ padded tiles, the volume the method was
   frozen on. With 200³ cores (240³ padded) the fixed K = 3600 candidate pool is diluted 1.73× and recall drops
   (`results/fullvolume/SMOKE_REPORT.md`).
4. **Face rejection** happens on the faces of the padded tile only (frozen v3 rejects candidates touching the
   faces of its input).
5. **Tile set:** GT-free, specimen-touching tiles only.
   - The specimen mask is built from the greyscale at step 8: Otsu, closing with a ball of radius 2, largest
     component, fill holes, dilation by 1.
   - A tile is selected when its core touches the mask: 277 of 810 tiles.
   - Post-hoc check: all 29,361 GT labels lie inside the selected cores.
6. **Queue order:** tiles whose padded coarse specimen fraction is ≥ 0.5 run first (163 tiles, grid order). Edge
   tiles run last, by decreasing fraction.
7. **Ownership:** an instance belongs to the tile whose core contains floor(centroid).
8. **Stitch:**
   - Tiles are merged sequentially in grid order.
   - Instances from different tiles are matched one-to-one at IoU ≥ 0.5, and the larger one is kept.
   - Contested voxels go to the tile whose core contains them.
   - **Orphan rescue is on:** a non-owned instance whose centroid lies in a selected core and whose voxels are
     still ≥ 50 % free is written into the free voxels only.
9. **Evaluation:** streaming evaluate_v2 (seed-ignore, dominant rule), valid GT ≥ 1000 voxels with its bbox inside
   the selected cores. Two GT sets:
   - `inside`: `labelsInside.txt` (9,278), the spec;
   - `listed`: `labelsInside` + `labelsOutside` (24,213).
10. **GT label counts (EFRGP01_00):**
    - 29,361 labels (every grain in the specimen) = 9,278 inside + 14,935 outside + 5,148 unlisted.
    - "Inside" is an inclined slab (the shear-band ROI): its xy position moves with z.
    - The unlisted grains are boundary grains: 4,300 lie more than 320 voxels from the axis (membrane shell) and
      1,160 are centred in the end caps; the two groups overlap.
