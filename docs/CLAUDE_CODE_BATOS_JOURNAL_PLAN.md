# BA-TOS Journal Extension — Experimental Plan for Claude Code

## Context

We are extending the SIBGRAPI 2026 paper:

**BA-TOS: Boundary-Aware Tree-of-Shapes Markers for 3D Grain Instance Segmentation in Micro-CT**

The goal is **not** to abandon BA-TOS or replace it with a completely different method.

The new objective is to make BA-TOS more mature and suitable for a **journal submission**, with:

1. stronger modern baselines;
2. foundation-model comparisons;
3. robustness and generalization experiments;
4. possible improvements to BA-TOS;
5. larger-scale/full-volume evaluation;
6. morphological usefulness of the segmented grains;
7. reproducible experiment scripts and logs.

The original SIBGRAPI BA-TOS implementation must be treated as a **frozen baseline**.

---

## Critical Rule

**DO NOT modify the original BA-TOS baseline code or overwrite its outputs.**

Create separate experiment folders and separate versions.

Suggested naming:

```text
BA-TOS_SIBGRAPI
BA-TOS_v2
BA-TOS_fullvolume
BA-TOS_hybrid
```

The published SIBGRAPI configuration must remain reproducible exactly as it was.

---

## Data

The original micro-CT volumes are:

```text
EFRGP01_00.tif
EFRGP01_01.tif
```

They are large TIFF volumes (~6 GB each).

On macOS, the intended local path is:

```text
~/Documents/BA_TOS_DATA/
```

Expected files:

```text
~/Documents/BA_TOS_DATA/EFRGP01_00.tif
~/Documents/BA_TOS_DATA/EFRGP01_01.tif
```

Do not run heavy experiments directly from the SMB/network share.

---

## First Task: Inspect the Existing BA-TOS Project

Before installing or changing anything:

1. locate the existing BA-TOS source code;
2. identify the exact scripts used for:
   - crop extraction;
   - preprocessing;
   - Tree of Shapes construction;
   - marker selection;
   - watershed;
   - evaluation;
   - metric computation;
3. identify the exact 12 crops used in the SIBGRAPI paper;
4. identify the frozen parameters;
5. identify the ground-truth files and their format;
6. identify the existing output format used for instance labels;
7. identify the Python environment / dependencies;
8. document all of the above in:

```text
docs/BA_TOS_BASELINE_REPRODUCTION.md
```

Do not refactor the baseline yet.

---

## Suggested Repository Structure

```text
BA_TOS_JOURNAL/
│
├── data/
│   ├── raw/
│   ├── crops/
│   ├── ground_truth/
│   └── metadata/
│
├── baselines/
│   ├── batos_sibgrapi/
│   ├── otsu_watershed/
│   ├── mintree_watershed/
│   ├── cellpose/
│   ├── nninteractive/
│   ├── sam_med3d/
│   ├── medical_sam3/
│   └── locate_anything/
│
├── experiments/
│   ├── batos_v2/
│   ├── batos_fullvolume/
│   └── batos_hybrid/
│
├── results/
│   ├── masks/
│   ├── metrics/
│   ├── runtime/
│   ├── memory/
│   └── figures/
│
├── scripts/
├── docs/
├── logs/
└── environments/
```

---

# Phase 1 — Reproduce BA-TOS SIBGRAPI

The first milestone is **exact reproduction** of the SIBGRAPI results.

Do not start foundation-model comparisons before this is working.

For every one of the 12 original crops, reproduce and save:

```text
prediction
TP
FP
FN
precision
recall
F1 @ IoU 0.25
F1 @ IoU 0.50
F1 @ IoU 0.75
number of predicted instances
number of GT instances
over-segmentation ratio
runtime
peak memory if possible
```

Also reproduce the reported aggregate values.

The frozen baseline should be stored under:

```text
results/batos_sibgrapi/
```

---

# Phase 2 — Cellpose First

Cellpose is the first modern baseline to test because it is relatively easy to run and supports 3D inference.

Primary models to test:

```text
cpsam_v2
cpdino
cpdino-vitb
```

## Important comparison rule

Use the **same 12 crops**, same ground truth, and same object-level evaluation protocol as BA-TOS.

Do not manually tune each crop.

Initial pass:

- default model settings;
- one fixed configuration per model;
- no per-crop optimization using GT.

After the default run, a small controlled sensitivity analysis may be performed.

Suggested variants:

```text
flow3D_smooth = 0, 1, 2
diameter = default and a small fixed set only if needed
```

Record:

```text
runtime
peak memory
number of instances
precision
recall
F1 @ 0.25
F1 @ 0.50
F1 @ 0.75
AP50
AP75
mean matched IoU
```

Save results under:

```text
results/cellpose/
```

---

# Phase 3 — Additional 3D Foundation / Promptable Models

After Cellpose works, evaluate selected 3D models.

Priority order:

1. nnInteractive
2. SAM-Med3D / SAM-Med3D-turbo
3. Medical-SAM3 or an equivalent 3D medical segmentation foundation model

These models are not necessarily directly comparable if they require prompts.

Therefore distinguish:

```text
automatic
oracle-prompt
BA-TOS-prompted
```

## Oracle-prompt experiment

Prompts derived from ground truth are allowed only as an **upper-bound experiment**.

Label them clearly:

```text
oracle prompt upper bound
```

Never compare an oracle-prompt result to automatic BA-TOS as if both used the same information.

## BA-TOS-prompted experiment

A key journal experiment is:

```text
micro-CT
  ↓
Tree of Shapes / BA-TOS
  ↓
automatic point/box/marker prompts
  ↓
3D foundation model
  ↓
refined grain instances
```

This could become a strong BA-TOS hybrid contribution.

---

# Phase 4 — LocateAnything

Project of interest:

```text
https://research.nvidia.com/labs/lpr/locate-anything/
```

Important:

**LocateAnything is not a native 3D instance segmentation model.**

Treat it as a 2D localization / grounding experiment.

Use it on selected 2D slices from the 3D volumes.

Suggested fixed slice positions:

```text
z = 50
z = 75
z = 100
z = 125
z = 150
```

Use fixed prompts, for example:

```text
grain
individual grain
granular particle
particle
```

Do not change the prompt per image after inspecting GT.

Measure:

```text
box recall
box precision
box IoU
missed grains
duplicate detections
```

Possible later hybrid:

```text
LocateAnything
  ↓
2D boxes / points
  ↓
SAM-like segmentation
  ↓
slice association / 3D reconstruction
  ↓
3D instances
```

This is exploratory and should not replace the core 3D comparison.

---

# Phase 5 — Improve BA-TOS

Only after the external baselines are measured should BA-TOS itself be modified.

## 1. Adaptive candidate capacity K

The SIBGRAPI paper showed sensitivity to insufficient candidate capacity.

Investigate whether `K` can become data-adaptive instead of fixed.

Possible signals:

```text
crop volume
foreground fraction
number of hierarchy nodes
candidate density
intensity statistics
estimated number of particles
```

Do not tune K independently for each test crop using ground truth.

## 2. Better support-collapse safeguard

Investigate a continuous or more general support-quality score.

Possible features:

```text
median support size
Q25 / Q75
support-size dispersion
duplicate-marker rate
number of retained nodes
foreground-to-marker ratio
hierarchical depth statistics
```

Any new gate must be frozen on development data before final evaluation.

## 3. Full-volume / tiled BA-TOS

A major journal extension should investigate larger regions or full volumes.

Potential strategy:

```text
large volume
  ↓
overlapping 3D tiles
  ↓
halo around each tile
  ↓
BA-TOS per tile
  ↓
retain instances by core ownership / centroid
  ↓
merge overlapping predictions
  ↓
one-to-one IoU reconciliation
```

Need to prevent boundary rejection from incorrectly deleting grains that only touch an artificial tile boundary.

## 4. Hybrid BA-TOS + foundation model

Potential contribution:

```text
BA-TOS markers / regions
  ↓
automatic prompts
  ↓
3D foundation model refinement
```

Research question:

> Can morphology-derived hierarchical markers provide effective automatic prompts for generalist 3D foundation models in granular micro-CT?

---

# Phase 6 — Journal-Grade Evaluation

All methods should use the same evaluation framework whenever methodologically valid.

Core segmentation metrics:

```text
Precision
Recall
F1 @ IoU 0.25
F1 @ IoU 0.50
F1 @ IoU 0.75
AP50
AP75
mean matched IoU
over-segmentation
under-segmentation
number of predicted instances
```

Robustness summaries:

```text
mean
standard deviation
median
minimum
95% bootstrap confidence interval
```

Compute both:

```text
pooled-object metrics
crop-macro metrics
```

Do not mix them in the same interpretation without labeling them.

---

# Morphological Accuracy

The journal paper should not evaluate only segmentation F1.

For matched grains, compare physical / morphological descriptors.

Suggested measurements:

```text
volume error
equivalent diameter error
centroid error
surface area error
sphericity error
elongation / aspect-ratio error
orientation error
```

This is important because the downstream goal is quantitative grain analysis.

---

# Runtime and Resource Evaluation

For every method record:

```text
wall-clock runtime
device
CPU / GPU / MPS
peak VRAM or unified memory where available
input size
model name
model version
software version
```

Create one machine-readable CSV/Parquet table for all runs.

Example columns:

```text
experiment_id
method
model
volume
crop
device
parameters
runtime_seconds
peak_memory_mb
tp
fp
fn
precision
recall
f1_025
f1_050
f1_075
ap50
ap75
mean_matched_iou
n_pred
n_gt
```

---

# Data Splitting / Generalization

The journal should not rely only on the same 12 development crops.

Preferred structure:

```text
Development:
original BA-TOS crops / known scans

Test:
new crops not used for parameter design

External:
another scan and preferably another material or acquisition condition
```

Parameters must be frozen before final test evaluation.

---

# Reproducibility Requirements

Every experiment must be reproducible from a command.

Create:

```text
scripts/run_batos_baseline.py
scripts/run_cellpose.py
scripts/run_evaluation.py
scripts/run_all_baselines.py
scripts/aggregate_results.py
scripts/morphology_analysis.py
```

Each run should save:

```text
command
timestamp
git commit
environment
parameters
input files
output files
metrics
```

Prefer YAML configuration files instead of hard-coded experiment settings.

Example:

```text
configs/
  batos_sibgrapi.yaml
  cellpose_cpsam_v2.yaml
  cellpose_cpdino.yaml
```

---

# Immediate Milestone on the Mac

Do only the following first:

1. verify the two TIFF files locally;
2. inspect the Mac hardware;
3. inspect Python/Conda;
4. locate and inspect the existing BA-TOS code;
5. reproduce one original BA-TOS crop if possible;
6. install Cellpose in a separate environment;
7. run one 3D crop through `cpsam_v2`;
8. confirm the output label volume;
9. run the exact BA-TOS evaluator on that output;
10. only then automate the 12-crop benchmark.

Do **not** try to process the full ~6 GB TIFF volume with Cellpose as the first experiment.

Start from the same 200³ crops used by BA-TOS.

---

# Deliverables Expected from Claude Code

For each step, Claude Code should:

1. inspect before modifying;
2. explain what file/script is being changed;
3. preserve the original BA-TOS implementation;
4. avoid deleting existing results;
5. use new experiment directories;
6. create reproducible commands;
7. validate array dimensions and label formats;
8. fail explicitly if GT or required metadata are missing;
9. never invent paths, metrics, or results;
10. maintain a running log in:

```text
docs/EXPERIMENT_LOG.md
```

---

# First Action

Start by inspecting the project and producing:

```text
docs/BA_TOS_BASELINE_REPRODUCTION.md
```

That document must identify:

- entry-point scripts;
- input data paths;
- crop definitions;
- GT paths;
- preprocessing;
- ToS implementation;
- frozen SIBGRAPI parameters;
- evaluation script;
- expected SIBGRAPI metrics;
- dependencies;
- exact command needed to reproduce one crop.

Do not modify the method until this reproduction document is complete.
