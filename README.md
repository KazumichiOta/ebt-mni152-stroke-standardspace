[![DOI (concept)](https://zenodo.org/badge/DOI/10.5281/zenodo.19472368.svg)](https://doi.org/10.5281/zenodo.19472368)

# EBT–MNI152 Stroke Standard-Space Comparison

# 1. Overview

This repository contains released derived data and reproducibility scripts related to our study,

*Standard-space selection in stroke magnetic resonance imaging: an elderly brain template reduces deformation bias and preserves lesion geometry*

This repository is intentionally separate from the main EBT repository, **`elderly-brain-template`**, which distributes the elderly brain template itself. The present repository focuses on the **stroke-specific comparison between EBT and MNI152** and its reproducibility package.

The study evaluates standard-space selection in lesion-bearing stroke MRI using 603 3T T1-weighted MRI scans from the ATLAS R2.0 Training Dataset. MNI152 and EBT normalization were performed under identical ANTs rigid→affine→SyN registration settings.

Specifically, this repository includes:

- released per-subject registration and whole-brain signed deformation metrics for **n = 603** included 3T ATLAS R2.0 cases,
- released lesion round-trip geometry metrics in native space,
- released lesion-wise log-Jacobian summaries and lesion volumes,
- a mask-definition sensitivity table based on subject-derived brain masks warped to each template space,
- subject-level selection and distance tables for the supplementary lentiform-anchored lesion-location analysis,
- the included subject list,
- scripts for regenerating key released derived data from a local ATLAS-based working directory,
- scripts for reproducing the main manuscript tables, figures, and key supplementary analyses from the released data.

ATLAS R2.0 subject-level MRI images and lesion masks are **not redistributed** in this repository.

---

# 2. Repository structure

```text
ebt-mni152-stroke-standardspace/
├─ analysis/
│  ├─ build_release_metrics_tsv.py
│  ├─ compute_lentiform_lesion_dist_native_roundtrip.py
│  ├─ compute_lesion_logjac_volume_T1_MNI_EBT.py
│  ├─ compute_ventricle_brain_fraction.py
│  ├─ make_all_from_release.R
│  ├─ recompute_lesion_roundtrip_dice_assd.py
│  └─ select_lentiform_subset_from_full_cohort.py
├─ data/
│  ├─ bg_CC_gradCC_subjectmask_warp_table.tsv
│  ├─ lentiform_lesion_dist_native_roundtrip.tsv
│  ├─ lentiform_selection_all_subjects.tsv
│  ├─ lesion_logjac_volume_T1_MNI_EBT.tsv
│  ├─ lesion_roundtrip_dice_assd_native.tsv
│  ├─ release_per_subject_metrics.tsv
│  └─ subjects_included.xlsx
└─ README.md
```

Generated figures and statistical outputs can be written to a local `outputs/` directory when running the reproducibility scripts. The generated `outputs/` directory does not need to be version-controlled.

---

# 3. Released data files

## 3.1 `data/release_per_subject_metrics.tsv`

**Description**

Released per-subject registration and whole-brain metrics for the **603 included 3T ATLAS R2.0 cases**.

**Main contents**

- lesion-excluded background correlation coefficient (CC) for MNI152 and EBT,
- gradient-CC for MNI152 and EBT,
- whole-brain signed log-Jacobian summaries,
- native-space lesion volume,
- ventricular fraction,
- scanner manufacturer,
- image-resolution metadata.

Whole-brain log-Jacobian summaries include mean, median, p5, and p95 values and characterize the signed distribution of contraction and expansion under each template condition.

The file is used to reproduce the principal whole-brain comparisons, ventricular-fraction analysis, lesion-volume supplementary regression analyses, and related figures.

---

## 3.2 `data/lesion_roundtrip_dice_assd_native.tsv`

**Description**

Released per-subject lesion round-trip geometry metrics in native space.

For each template condition, the native-space lesion mask was transformed to template space and then mapped back to native space using the corresponding template-specific forward and inverse transformations.

**Main contents**

- Dice coefficient for the MNI152 condition,
- Dice coefficient for the EBT condition,
- average symmetric surface distance (ASSD; mm) for the MNI152 condition,
- ASSD for the EBT condition.

These data are used for:

- paired between-template comparisons,
- prespecified TOST equivalence analyses,
- round-trip lesion geometry figures,
- supplementary analyses relating round-trip boundary error to intralesional deformation.

---

## 3.3 `data/lesion_logjac_volume_T1_MNI_EBT.tsv`

**Description**

Released lesion-wise signed log-Jacobian summaries and lesion-volume measurements.

**Main contents**

- native-space lesion volume,
- template-space lesion volumes,
- intralesional mean log-Jacobian,
- intralesional median log-Jacobian,
- intralesional p5 log-Jacobian,
- intralesional p95 log-Jacobian,
- corresponding values under the MNI152 and EBT conditions.

The intralesional mean log-Jacobian values are used in the main manuscript and in supplementary analyses examining their relationship with round-trip ASSD and lesion volume.

---

## 3.4 `data/bg_CC_gradCC_subjectmask_warp_table.tsv`

**Description**

Supplementary mask-definition sensitivity data based on subject-derived brain masks.

Native-space HD-BET brain masks were warped to each template space and used to define an alternative lesion-excluded background region.

**Main contents**

- subject-mask background CC for MNI152,
- subject-mask background CC for EBT,
- subject-mask gradient-CC for MNI152,
- subject-mask gradient-CC for EBT.

This file is used to reproduce Supplementary Table S4.

---

## 3.5 `data/lentiform_selection_all_subjects.tsv`

**Description**

Subject-level audit table for reconstructing the supplementary lentiform-anchored lesion-location analysis from the full study cohort.

The selection criteria were:

1. a single contiguous native-space lesion, and
2. no lesion voxels within the bilateral SynthSeg-derived lentiform nucleus.

The lentiform nucleus was defined as the **bilateral putamen and pallidum** using SynthSeg labels:

- 12: left putamen,
- 13: left pallidum,
- 51: right putamen,
- 52: right pallidum.

A single contiguous lesion was defined as one 3D connected component using **26-connectivity**.

For the manuscript dataset, the reproduced selection flow was:

```text
Full cohort:                    603
Single contiguous lesion:      279
No lentiform lesion:           114
Final analysis subset:         114
```

The CST-overlap criterion was **not** used for this supplementary analysis.

The table provides the subject-level inclusion and exclusion status used to identify the final analysis subset.

For Supplementary Table S8, subjects with:

```text
final_include == TRUE
```

are retained.

---

## 3.6 `data/lentiform_lesion_dist_native_roundtrip.tsv`

**Description**

Released subject-level measurements for the supplementary lentiform-anchored lesion-location analysis.

The file contains the 3D Euclidean distance between the center of mass (COM) of the native-space bilateral lentiform nucleus and the lesion COM for:

- the original native-space lesion,
- the MNI152 round-tripped lesion,
- the EBT round-tripped lesion.

The lentiform nucleus is defined as the bilateral putamen and pallidum using SynthSeg labels:

```text
12, 13, 51, 52
```

Round-trip distance errors are calculated as:

```text
err_MNI = |dist_MNIround − dist_native|
err_EBT = |dist_EBTround − dist_native|
```

Supplementary Table S8 is reconstructed by combining this file with:

```text
data/lentiform_selection_all_subjects.tsv
```

and retaining subjects with:

```text
final_include == TRUE
```

The resulting final analysis subset contains **n = 114** subjects.

---

## 3.7 `data/subjects_included.xlsx`

**Description**

List of the **603 ATLAS R2.0 3T subjects** included in the final analysis cohort.

This file is provided for:

- cohort auditing,
- subject-level cross-checking,
- confirmation of the final included dataset.

It is not required for reproduction of the statistical analyses from the released TSV files.

---

# 4. Analysis scripts

## 4.1 `analysis/build_release_metrics_tsv.py`

**Purpose**

Generates the released per-subject registration and whole-brain metrics table from a local ATLAS-based working directory.

**Main output**

```text
release_per_subject_metrics.tsv
```

The script combines previously derived subject-level registration, deformation, lesion-volume, ventricular-fraction, scanner, and image-resolution variables into the public release table.

---

## 4.2 `analysis/recompute_lesion_roundtrip_dice_assd.py`

**Purpose**

Recomputes lesion round-trip Dice and ASSD in native space using template-specific forward and inverse transforms.

For each template condition:

```text
native lesion
    ↓
template space
    ↓
native space
```

The original native-space lesion mask is then compared with the round-tripped lesion mask.

**Main output**

```text
lesion_roundtrip_dice_assd_native.tsv
```

---

## 4.3 `analysis/compute_lesion_logjac_volume_T1_MNI_EBT.py`

**Purpose**

Computes lesion-wise signed log-Jacobian summaries and lesion volumes for the MNI152 and EBT conditions.

**Main output**

```text
lesion_logjac_volume_T1_MNI_EBT.tsv
```

The rigid component is excluded from log-Jacobian calculations because rigid transformations are volume-preserving.

---

## 4.4 `analysis/compute_ventricle_brain_fraction.py`

**Purpose**

Computes ventricular fraction from native-space SynthSeg segmentation.

Ventricular volume is defined using SynthSeg labels corresponding to:

- bilateral lateral ventricles,
- bilateral inferior lateral ventricles,
- third ventricle,
- fourth ventricle.

The labels used are:

```text
4, 5, 14, 15, 43, 44
```

Intracranial volume (ICV) is defined as the total volume of voxels with a nonzero SynthSeg label.

Ventricular fraction is calculated as:

```text
ventricular fraction = ventricular volume / ICV
```

The resulting ventricular fraction is used as an imaging-based proxy of atrophy and ventricular enlargement.

---

## 4.5 `analysis/select_lentiform_subset_from_full_cohort.py`

**Purpose**

Reconstructs the subject subset used for the supplementary lentiform-anchored lesion-location analysis directly from the full study cohort.

The script applies the actual selection criteria used for this analysis:

1. single contiguous native-space lesion,
2. no lesion voxels within the bilateral SynthSeg-derived lentiform nucleus.

The single-lesion criterion uses 3D connected-component labeling with **26-connectivity**.

For the manuscript dataset, the script reproduces:

```text
603 → 279 → 114
```

corresponding to:

```text
603 total cases
279 cases with a single contiguous lesion
114 cases with no lesion voxels within the lentiform nucleus
```

**Main output**

```text
lentiform_selection_all_subjects.tsv
```

This output provides an auditable subject-level record of inclusion and exclusion across the full cohort.

The final analysis subset is identified using:

```text
final_include == TRUE
```

---

## 4.6 `analysis/compute_lentiform_lesion_dist_native_roundtrip.py`

**Purpose**

Computes the 3D Euclidean distance between the COM of the native-space bilateral lentiform nucleus and the lesion COM.

Distances are calculated for:

- the original native-space lesion,
- the MNI152 round-tripped lesion,
- the EBT round-tripped lesion.

The lentiform nucleus is defined as the bilateral putamen and pallidum using SynthSeg labels:

```text
12, 13, 51, 52
```

**Main output**

```text
lentiform_lesion_dist_native_roundtrip.tsv
```

Round-trip distance errors are calculated as:

```text
err_MNI = |dist_MNIround − dist_native|
err_EBT = |dist_EBTround − dist_native|
```

For Supplementary Table S8, these measurements are restricted to the single-lesion cases without lesion involvement of the lentiform nucleus.

---

## 4.7 `analysis/make_all_from_release.R`

**Purpose**

Reproduces the principal manuscript tables, figures, and key supplementary statistical analyses directly from the released TSV files in `data/`.

The script uses the following primary released inputs:

```text
release_per_subject_metrics.tsv
lesion_roundtrip_dice_assd_native.tsv
lesion_logjac_volume_T1_MNI_EBT.tsv
```

and, when available:

```text
bg_CC_gradCC_subjectmask_warp_table.tsv
```

for the mask-definition sensitivity analysis.

For the supplementary lentiform-anchored lesion-location analysis, the script uses:

```text
lentiform_selection_all_subjects.tsv
lentiform_lesion_dist_native_roundtrip.tsv
```

The two files are combined at the subject level, and subjects satisfying:

```text
final_include == TRUE
```

are retained to reconstruct the final **n = 114** analysis subset and Supplementary Table S8.

**Typical outputs include**

```text
Table2.tsv
Table3.tsv
Fig2A_CC_vs_lesionVol.pdf
Fig2B_wholebrain_mean_logJ.pdf
Fig3A_dDice_margins.pdf
Fig3B_dASSD_margins.pdf
Fig3C_lesion_mean_logJ.pdf
Fig4_dCC_vs_ventricle_fraction.pdf
S4_subjectmask_*.tsv
S5_lesionVol_regression.tsv
S6_absMeanLogJ_vs_ASSD.tsv
S7_intralesional_dlogJ_tertiles.tsv
S8_lentiform_lesion_distance.tsv
S9_ventricle_fraction_dCC.tsv
```

The exact output filenames may vary slightly with the release version of the script.

---

# 5. Reproducing the main manuscript analyses from released data

From the repository root, run:

```bash
Rscript analysis/make_all_from_release.R \
  --root . \
  --outdir outputs
```

This command reproduces the principal statistical analyses directly from the released derived data without requiring access to the original ATLAS R2.0 MRI files.

The resulting analyses include:

- MNI152 versus EBT background-CC comparison,
- gradient-CC comparison,
- whole-brain signed log-Jacobian comparisons,
- lesion round-trip Dice comparison,
- lesion round-trip ASSD comparison,
- TOST equivalence testing,
- intralesional mean log-Jacobian comparison,
- ventricular-fraction versus ΔCC regression,
- lesion-volume supplementary regression models,
- association between the absolute value of intralesional mean logJ and round-trip ASSD,
- intralesional logJ analyses by lesion-volume tertile,
- mask-definition sensitivity analysis,
- reconstruction of the supplementary lentiform-anchored lesion-location analysis.

Supplementary Table S8 is reconstructed from:

```text
data/lentiform_selection_all_subjects.tsv
data/lentiform_lesion_dist_native_roundtrip.tsv
```

by retaining subjects with:

```text
final_include == TRUE
```

The expected final subset is:

```text
n = 114
```

---

# 6. Regenerating derived data from a local ATLAS working directory

The following scripts are provided for transparency and require access to a locally prepared ATLAS R2.0-derived working directory and the corresponding preprocessing and registration outputs.

The repository does not contain those source MRI files.

---

## 6.1 Whole-brain released metrics

```bash
python analysis/build_release_metrics_tsv.py \
  --root /path/to/ATLAS_workdir
```

**Output**

```text
release_per_subject_metrics.tsv
```

---

## 6.2 Lesion round-trip Dice and ASSD

```bash
python analysis/recompute_lesion_roundtrip_dice_assd.py \
  --root /path/to/ATLAS_workdir
```

**Output**

```text
lesion_roundtrip_dice_assd_native.tsv
```

---

## 6.3 Lesion-wise log-Jacobian summaries and lesion volumes

```bash
python analysis/compute_lesion_logjac_volume_T1_MNI_EBT.py \
  --root /path/to/ATLAS_workdir \
  --mni-t1 /path/to/MNI152_T1_1mm_brain.nii.gz \
  --ebt-t1 /path/to/Elderly_brain_T1_template_1mm.nii.gz
```

**Output**

```text
lesion_logjac_volume_T1_MNI_EBT.tsv
```

---

## 6.4 Ventricular fraction

```bash
python analysis/compute_ventricle_brain_fraction.py \
  --root /path/to/ATLAS_workdir
```

Ventricular volume is defined from SynthSeg labels:

```text
4, 5, 14, 15, 43, 44
```

corresponding to the bilateral lateral ventricles, bilateral inferior lateral ventricles, third ventricle, and fourth ventricle.

Intracranial volume (ICV) is defined as the total volume of all voxels with a nonzero SynthSeg label.

The resulting ventricular fraction can be incorporated into the released per-subject metrics table.

---

## 6.5 Lentiform-analysis subject selection

```bash
python analysis/select_lentiform_subset_from_full_cohort.py \
  --root /path/to/ATLAS_workdir
```

For the manuscript dataset, the expected selection flow is:

```text
Full cohort                          : 603
Single contiguous lesion            : 279
+ no lesion within lentiform nucleus: 114
Final eligible subset               : 114
```

The lentiform nucleus is defined as the bilateral putamen and pallidum using SynthSeg labels:

```text
12, 13, 51, 52
```

A single contiguous lesion is defined using 3D connected-component labeling with 26-connectivity.

A successful reproduction should end with:

```text
[OK] Selection reproduced: 603 -> 279 -> 114.
```

**Output**

```text
lentiform_selection_all_subjects.tsv
```

---

## 6.6 Lentiform-to-lesion round-trip distance

```bash
python analysis/compute_lentiform_lesion_dist_native_roundtrip.py \
  --root /path/to/ATLAS_workdir
```

The script computes subject-level native, MNI152 round-trip, and EBT round-trip lentiform-to-lesion COM distances.

**Output**

```text
lentiform_lesion_dist_native_roundtrip.tsv
```

For Supplementary Table S8, these measurements are combined with:

```text
lentiform_selection_all_subjects.tsv
```

and restricted to subjects with:

```text
final_include == TRUE
```

The expected final analysis subset is **n = 114**.

---

# 7. Interpretation of the lentiform supplementary analysis

The supplementary lentiform analysis was designed as an additional check of relative lesion location after round-trip transformation.

It does **not** test CST overlap and does not use CST overlap as an inclusion criterion.

The analysis is restricted to cases with:

1. a single contiguous lesion, and
2. no lesion voxels within the bilateral lentiform nucleus.

The purpose of the analysis is to assess whether switching standard spaces produces gross changes in the spatial relationship between a lesion and a stable deep gray-matter landmark.

The lentiform nucleus was therefore used as an anatomical anchor, and the native-space lentiform-to-lesion COM distance was compared with the corresponding distance after MNI152 and EBT round-trip transformation.

The final analysis subset contained **114 cases**.

For the manuscript dataset, the reconstructed group-level values are approximately:

```text
Native distance             ≈ 51.15 ± 19.26 mm
MNI152 round-trip distance  ≈ 51.16 ± 19.25 mm
EBT round-trip distance     ≈ 51.15 ± 19.26 mm
MNI152 distance error       ≈ 0.04 ± 0.05 mm
EBT distance error          ≈ 0.04 ± 0.05 mm
Δerror (EBT − MNI152)       ≈ 0.0007 mm
95% CI                      ≈ −0.0089 to 0.0103 mm
paired t-test p             ≈ 0.884
```

These measurements are provided as an additional check of preservation of relative lesion location under the two standard-space conditions.

---

# 8. Data redistribution and source-data scope

This repository does **not** redistribute ATLAS R2.0 subject-level MRI images or lesion masks.

Only the following are provided:

- derived numerical metrics,
- subject-level derived analysis tables,
- supporting subject identifiers,
- reproducibility scripts.

Released files in `data/` are derived from the **ATLAS R2.0 Training Dataset** and are provided for transparency and reproducibility of the present study.

The included subject list documents the final study cohort but does not replace access to the original ATLAS R2.0 dataset.

Users who wish to regenerate released derived data from source images must obtain ATLAS R2.0 independently and comply with its original terms of use.

ATLAS R2.0 should be cited as:

Liew SL, Tavenner BP, Donnelly MR, Zavaliangos-Petropulu A, Jeong JN, Barisano G, et al.  
*A large, curated, open-source stroke neuroimaging dataset to improve lesion segmentation algorithms.*  
**Scientific Data.** 2022;9:320.  
https://doi.org/10.1038/s41597-022-01401-7

---

# 9. Software environment

The released scripts were developed using:

- Python 3.11
- R 4.4
- NumPy
- SciPy
- pandas
- NiBabel
- ANTs
- ANTsPy
- FSL
- HD-BET
- SynthSeg / FreeSurfer

Please refer to Supplementary Table S1 of the associated manuscript for the full software list and exact versions used in the study.

---

# 10. License

## 10.1 Code

All code in `analysis/` is released under the **MIT License**.

### MIT License

Copyright (c) 2026 Kazumichi Ota

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## 10.2 Released derived data

Released derived data in `data/` are provided under **CC BY 4.0**, unless otherwise noted.

You are free to:

- share,
- adapt,

under the terms of the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license, provided that appropriate credit is given.

---

## 10.3 Third-party source data

ATLAS R2.0 subject-level MRI images and lesion masks are **not** redistributed in this repository.

Use of the source dataset remains subject to the original ATLAS R2.0 terms.

---

# 11. Citation

If you use these released data or scripts, please cite the associated article and repository DOI.

**Article (provisional; update when finalized)**

Ota K, Nakazato Y, Oyama G.  
*Standard-space selection in stroke magnetic resonance imaging: an elderly brain template reduces deformation bias and preserves lesion geometry.*

**Repository DOI (this version)**

https://doi.org/<ZENODO_V1.1_DOI>

**Repository DOI (concept DOI; always latest)**

https://doi.org/10.5281/zenodo.19472368

---

# 12. Related repository

The elderly brain template itself is distributed separately through the main EBT repository:

**`elderly-brain-template`**

The present repository does not redistribute the EBT template and focuses specifically on the comparison between EBT and MNI152 in lesion-bearing stroke MRI.

---

# 13. Contact

**Corresponding author:** Kazumichi Ota  
**Email:** kota24@saitama-med.ac.jp

For questions, feedback, or reproducibility issues related to the released data or scripts, please contact the corresponding author.