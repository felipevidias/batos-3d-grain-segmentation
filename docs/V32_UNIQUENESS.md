# BA-TOS v3.2 (per-grain marker uniqueness) — design record, Phase 5h

Written 2026-11-09 **before any v3.2 run**. Design data: the 12 EFR dev crops (f = 1.0) and the synthetic scale set
(f = 0.7 / 1.4 / 2.0, `data/metadata/crops_synth_scale.json`). Test set 1 and test set 2 are not read during design.
Work in the git worktree `..\batos-v32`, branch `journal-extension-v32`.

## Motivation

- v3 on the EO interior crops of test set 1 has 0.35 markers per grain, and about 1/3 of the grains are unseeded.
  The cause is capacity: MIN candidates sit at the 1800/s³ cap on every crop.
- v3.1 (count constants floored at v2: K, tos_seed nodes, legacy markers) raised markers per grain to 2.62 at f = 2.0.
  It then over-segmented into sub-filter fragments (`results/v31_dev/REPORT.md`).
- So capacity alone cannot fix it. v3.2 = v3.1 + a uniqueness rule applied to the e11 ToS-MIN marker list.
  Everything downstream is the unchanged v3 tail: EDT-Otsu control matched to the new count, fixed-mask watershed,
  E23c gate, selective watershed.

## Hypotheses (`src/batos_v32/uniqueness.py`; s = r_granulo / 10)

| id | rule | new constants |
|---|---|---|
| H1 | One marker per basin of the EDT watershed of the Otsu mask. Mask = the e11 control/watershed Otsu mask (min component 500 s³, closing round(s)). Basin seeds = the EDT-Otsu control's own markers before count matching (EDT local maxima, greedy NMS 6 s). Per basin, keep the highest-priority ToS marker (ties broken by source row), as e11 does per internal component. Markers outside the mask are kept. | none (all reused from v3) |
| H2 | Greedy NMS on the final marker list in priority order, radius 0.8 r_granulo. | 0.8 (user-specified) |
| H3 | Merge markers whose exclusive ToS supports touch (6-adjacency; supports = the e10 nearest-selected-ancestor leaf labels, `batos_v32/supports.py`) and whose combined support is ≤ V_g = 4/3 π r_granulo³. Agglomerative: smallest combined pairs first, with union-find. Each group keeps its highest-priority marker. | V_g (user-specified) |

H1's basin seeding is not given in the user's specification. The choice above reuses the pipeline's own EDT-Otsu
marker rule, so no new constant is introduced. It is fixed here, before any result.

Mechanics (`src/batos_v32/pipeline.py`):
- The v3.1 prefix (scale → … → e11 ToS-MIN markers) runs once per crop.
- Then, per hypothesis: uniqueness rule, then tail.
- Hypothesis `none` must reproduce v3 (EFRGP01_00_crop04, f = 1) and v3.1 (EFRGP01_00_crop02_x2, f = 2)
  bit-for-bit (`identity_check.json`).

## Order and acceptance (fixed before running)

The hypotheses are tested one at a time, in the order H1 → H2 → H3. The first one that holds is frozen, and the
remaining ones are not run.

Primary metric: crop-macro F1@0.50 for BA-TOS_sel, seed-ignore, with no size filter. This is the same metric as in
Phase 5e/5f. A hypothesis **holds** if all of the following are true:

1. f = 1.0: crop-macro F1@0.50 ≥ 0.985.
2. f = 1.4 and f = 2.0: crop-macro F1@0.50 ≥ v3 − 0.005 on each (v3: 0.995 / 0.938). This is "holds or improves",
   with a 0.005 tolerance.
3. f = 0.7, a supplementary check because the uniqueness rule also acts at s ≤ 1: crop-macro F1@0.50 ≥ v3 − 0.005
   (v3: 0.975).

Reported per factor:
- markers per grain (all GT, the 5f definition, and per valid grain);
- fraction of valid grains seeded (≥ 1 marker inside) and split (≥ 2);
- F1@0.50 and F1@0.75, without a filter and with the own filter;
- pooled F1.

GT is used only for evaluation and diagnosis, never inside the method.

## If a hypothesis holds

1. Freeze it as `batos-paper-v3.2`: config, tag, and the SHA-256 of all executed code.
2. Cut test set 2: seed 20261109, the same GT-free rule as test set 1, 52 crops, no overlap with the dev crops or
   test set 1.
3. Run v3 and v3.2 on test set 2 once, and report per specimen.

If none holds: report and stop.
