# Dataset inventory — Pinzón et al. (Zenodo 8014905), six specimens

Written 2026-11-06 (Phase 1f) from headers and archive listings only. Nothing was extracted or loaded in full.
Root: `C:\Users\USER\Downloads\PhDUpload\PhDUpload\<SPEC>\`. Every specimen folder has the same subfolders:
`Greyscale`, `Label`, `Tracking`, `Particles`, `Contacts`, `Strain`, `MacroscopicStressStrain`, `ROI`.
The meaning of the specimen codes (F / M / O and R / S) is not stated in these files. It is not guessed here; see
Pinzón et al., *Open Geomechanics* 2025.

## Volumes and timesteps

| specimen | greyscale present | greyscale shape (z, y, x), dtype | label archive (7z, non-solid, LZMA2) | label members | timesteps |
|---|---|---|---|---|---|
| EFR | `EFRGP01_00.tif`, `_01.tif` | 1600 × 1400 × 1400, uint16, 1600 pages, 6.27 GB each | `EFR_labels.7z`, 2.62 GB packed / 345 GB unpacked | `EFRGP01_00.tif` … `EFRGP01_54.tif` (55, contiguous, ~6.27 GB each) | 55 (00–54) |
| EFS | `EFSGP01_00.tif`, `_01.tif` | same | `EFS_labels.7z`, 2.66 GB | `EFSGP01_00` … `_54` (55) | 55 |
| EMR | `EMRGP01_00.tif`, `_01.tif` | same | `EMR_labels.7z`, 1.32 GB | `EMRGP01_00` … `_30` (31) | 31 (00–30) |
| EMS | `EMSGP01_00.tif`, `_01.tif` | same | `EMS_labels.7z`, 1.30 GB | `EMSGP01_00` … `_30` (31) | 31 |
| EOR | `EORGP01_00.tif`, `_01.tif` | same | `EOR_labels.7z`, 1.26 GB | `EORGP01_00` … `_30` (31) | 31 |
| EOS | `EOSGP01_00.tif`, `_01.tif` | same | `EOS_labels.7z`, 1.26 GB | `EOSGP01_00` … `_30` (31) | 31 |

**Greyscale.** Only the timesteps `_00` and `_01` are present locally, for all six specimens.

**Labels.** Each archive holds one full label volume per timestep, with the same shape as the greyscale. The archive
is non-solid, so `_00` and `_01` can be extracted alone. EFR `_00` / `_01` were already extracted to
`C:\BA_TOS_DATA\gt_full\` (Phase 1b).

**Specimen geometry** (greyscale only, from `scripts/define_test_crops.py`): in all 12 scans the specimen is a
vertical cylinder with equivalent radius 350–352 voxels and its axis at about (x, y) = (707, 711). The constant-section
part spans z ≈ 104–1476; the end caps lie outside it.

**Grain counts** (number of rows in the first tracking file): EFR 29,362; EFS 29,424; EMR 19,980; EMS 19,972;
EOR 19,076; EOS 19,248. The M and O specimens have about one third fewer grains in the same specimen volume, so their
grains are larger than in the F specimens.

## Per-specimen text / TIFF products (EFR shown; the same layout everywhere)

| folder | files | content (first lines / header) |
|---|---|---|
| `Tracking\` | `<S>GP01_tt_tt+1.txt`, one per increment (EFR / EFS 54, others 30) | Tab-separated, header `Label Zpos Ypos Xpos Zdisp Ydisp Xdisp Fzz Fzy Fzx Fyz Fyy Fyx Fxz Fxy Fxx PSCC error iterations returnStatus deltaPhiNorm LabelDilate`. One row per label (row 0 = background); per-grain position, displacement and deformation gradient F from DVC tracking (spam-style). EFR 00→01: 29,363 rows. |
| `Particles\COM\` | `<S>GP01_tt_tt+1.txt` (54 / 30) | `# Z - Y - X`; per-label centre of mass (row 0 = background). |
| `Particles\Orientation\` | `<S>GP01_tt_Particles.txt` (55 / 31) | `# Z[deg] - Y[deg] - X[deg]`; per-grain orientation vector. |
| `Particles\Rotation\` | `<S>GP01_tt_tt+1_particleRotation.txt` | `# Mag [deg] - Orientation Z - Y - X`; per-grain incremental rotation. |
| `Contacts\ContactPoint\` | `<S>GP01_tt.txt` (55 / 31) | `# Grain A - Grain B - Z - Y - X`; contact pairs and contact-point position. |
| `Contacts\Orientation\` | `<S>GP01_tt.txt` | `# Grain A - Grain B - Z - Y - X - intervox`; contact normal and contact size in voxels. |
| `Contacts\Dynamics\` | `<S>GP01_tt_tt+1_{Created,…}.txt` (162 in EFR) | `# Grain A - Grain B - Z - Y - X - intervox`; contacts created / lost between increments. |
| `Strain\IncrementalFields(_Rotated)\` | `<S>GP01-tt-tt+1-{dev,vol}Strain.tif` | Coarse float32 strain fields (EFR 44 × 38 × 38, EMR 38 × 33 × 33). |
| `Strain\MacroStrain\` | `<S>GP01_macroStrain.txt` | `# Scan AxialStrain [%] VolumetricStrain [%]`, one row per scan. |
| `MacroscopicStressStrain\` | `<S>GP01_macroStressStrain.txt` (EFR 17,166 rows) | `# AxialStrain [%] - AxialStress [kPa] - DevStress [kPa] - MeanStress [kPa] - q/p [-]`. The first row is about 50 kPa mean stress for every specimen (confinement). |
| `ROI\` | `<S>GP01_labelsInside.txt`, `<S>GP01_labelsOutside.txt` | `# Labels`, then one label ID per line. EFR: 9,278 inside / 14,935 outside (EFS 8,834 / 15,511; EMR 6,060 / 9,269; EMS 6,197 / 9,855; EOR 5,243 / 10,093; EOS 6,791 / 8,690). A partition of the grain labels by the authors' region of interest, presumably a central sub-volume; the files do not define it. **Not used** for test-crop selection: the rule uses the greyscale membrane mask only. |

## Use in this project
- **Development:** 12 EFR crops (`crops_sibgrapi.json`, `crops_manifest.csv`).
- **Test (Phase 1f):** 52 crops, `crops_test.json` → `crops_manifest_test.csv`, cut to `C:\BA_TOS_DATA\crops_test\`:
  - EFR internal: 6 crops per scan;
  - EFS, EMR, EMS, EOR, EOS external: 4 crops per scan.
- **Later timesteps** (e.g. 27, 54) need greyscale volumes that are not present locally (`scripts/cut_timestep_crops.py`).
