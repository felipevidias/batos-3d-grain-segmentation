# BA-TOS: Boundary-Aware Tree-of-Shapes Markers for 3D Grain Instance Segmentation in Micro-CT

This repository contains research code, scripts, figures, and preliminary experimental material for a morphology-based pipeline for **3D grain instance segmentation in micro-computed tomography (micro-CT)** images.

The current approach uses **Tree-of-Shapes-based connected-operator processing** as a marker generation and refinement step before marker-controlled distance watershed segmentation.

> **Research status:** this repository is an active research draft. The method name, paper framing, and final experiments may still change after advisor feedback.

## Overview

Instance segmentation of individual grains in 3D micro-CT images is challenging because granular materials often contain densely packed particles, touching grains, dark inter-grain/background regions, boundary-connected structures, local intensity variation, and unstable watershed markers.

The proposed pipeline aims to improve the marker generation step used by watershed segmentation. Instead of relying only on local distance-transform maxima or Min-tree cores, the method uses Tree-of-Shapes-based processing to extract candidate internal regions, rejects dark regions connected to the crop boundary, and keeps refined internal candidates as 3D markers.

## Current pipeline

1. **Conservative 3D flattening** of the input crop;
2. **Tree-of-Shapes-based candidate extraction**;
3. **Rejection of dark components connected to the crop boundary**;
4. **Refinement of internal candidates into 3D markers**;
5. **Marker-controlled distance watershed**;
6. **Conservative size-aware post-filtering** of small predicted instances.

The Tree of Shapes is not used here as the final segmentation method. It is used as a **topology-guided marker generation/refinement mechanism** before watershed segmentation.

## Compared methods

The current experiments compare:

1. **Lysandre Min-tree + seeded watershed** — baseline inspired by the previous Min-tree core extraction + seeded watershed solution explored in Lysandre Macke's work on 3D CT grain segmentation.
2. **Otsu + distance watershed** — classical baseline based on global Otsu thresholding, distance transform, and watershed.
3. **ToS-based marker refinement + watershed** — proposed pipeline using boundary-aware Tree-of-Shapes-based internal marker refinement before distance watershed.

## Preliminary results

The current evaluation was performed on two consecutive scans, `EFRGP01_00` and `EFRGP01_01`. For each scan, six representative 3D crops of size `200 x 200 x 200` voxels were used, resulting in **12 crops**. The same parameter configuration was kept for both scans.

| Method | P@0.50 | R@0.50 | F1@0.50 | FP | FN |
|---|---:|---:|---:|---:|---:|
| Lysandre Min-tree + watershed | 0.560 | 0.386 | 0.436 | 509 | 1914 |
| Otsu + distance watershed | 0.840 | 0.961 | 0.895 | 607 | 115 |
| ToS-based marker refinement + watershed | 0.980 | 0.967 | 0.974 | 39 | 77 |

These preliminary results suggest that the Min-tree + watershed baseline tends to under-recover valid particles, Otsu + watershed obtains high recall but produces many false positives and fragmented labels, and the ToS-based marker refinement pipeline preserves high recall while reducing false positives.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── requirements.txt
├── src/batos/
├── scripts/
├── sample_data/
├── paper/
├── results/
└── docs/
```

## Sample data

The folder `sample_data/` contains two **small synthetic 3D crops** that can be committed to GitHub and used only to test the repository workflow.

Important:

- these demo crops are **not** Gustavo's real data;
- they are **not** used for the scientific results reported in the paper draft;
- they exist only to test loading, preprocessing, metrics, and scripts without downloading large micro-CT volumes.

## Data availability

The raw micro-CT volumes and reference labels are **not included** in this repository because of file size and distribution constraints. Users should download the original data from the official data source and adjust paths in the scripts.

Large raw data files such as `.tif`, `.tiff`, `.7z`, `.zip`, and full intermediate result folders should not be committed to this repository.

## Installation

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Main dependencies:

```text
numpy
scipy
scikit-image
tifffile
pillow
matplotlib
pandas
py7zr
higra
```

## Basic usage

Count labels in a test crop:

```bash
python src/batos/count_tiff_labels.py sample_data/synthetic_crop01_label.tif
```

Evaluate a label volume against itself as a sanity check:

```bash
python src/batos/batos_instance_metrics.py \
  --gt sample_data/synthetic_crop01_label.tif \
  --pred sample_data/synthetic_crop01_label.tif \
  --out-dir results/demo_eval_instance \
  --min-gt-area 50 \
  --thresholds 0.50 \
  --name synthetic_self_check
```

## Figures

Experimental visual comparisons should use real dataset images. AI-generated or schematic figures should only be used to explain concepts, not as experimental evidence.

## Paper draft

The current paper framing is:

> Boundary-aware Tree-of-Shapes-based marker refinement for 3D grain instance segmentation in micro-CT.

For double-blind submission, remove author names, affiliations, e-mails, acknowledgments, identifying metadata, and self-identifying repository links.

## Current limitations

- Evaluation is currently based on representative 3D crops, not full-volume segmentation.
- Experiments are focused on two consecutive scans from the same EFRGP01 sequence.
- Full-volume or slab-based processing still needs to be implemented.
- Label merging across overlapping 3D tiles remains future work.
- The pipeline and scripts are still being cleaned for reproducibility.

## Planned next steps

- Organize scripts into a cleaner command-line interface.
- Document exact data preparation steps.
- Add full-volume or overlapping-slab processing.
- Implement label reconciliation across tile borders.
- Evaluate additional scans and material conditions.
- Refine paper figures and captions.
- Prepare a clean reproducibility package after advisor feedback.

## Citation

A formal citation will be added if the associated paper is submitted or accepted. For now, please cite this repository as an ongoing research draft.

## Authors

- Felipe Vilhena Dias — IMScience Lab, PUC Minas
- Silvio Jamil F. Guimarães — IMScience Lab, PUC Minas
- Yukiko Kenmochi — GREYC Laboratory, University of Caen Normandy

## License

License pending. This repository is currently a private research draft. A final license will be selected after checking third-party code dependencies, dataset constraints, and advisor approval.
