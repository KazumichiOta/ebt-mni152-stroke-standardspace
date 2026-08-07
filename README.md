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
- a released summary table for the supplementary lentiform-to-lesion round-trip analysis,
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
│  ├─ compute_lesion_logjac_volume_T1_MNI_EBT.py
│  ├─ recompute_lesion_roundtrip_dice_assd.py
│  ├─ compute_ventricle_brain_fraction.py
│  ├─ compute_lentiform_lesion_dist_native_roundtrip.py
│  ├─ select_lentiform_subset_from_full_cohort.py
│  └─ make_all_from_release.R
├─ data/
│  ├─ bg_CC_gradCC_subjectmask_warp_table.tsv
│  ├─ lentiform_lesion_dist_noLentLesion_table.tsv
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

## 3.5 `data/lentiform_lesion_dist_noLentLesion_table.tsv`

**Description**

Released summary table for the supplementary lentiform-to-lesion center-of-mass (COM) round-trip analysis.

The underlying analysis was restricted to subjects with:

1. a single contiguous lesion, and
2. no lesion voxels within the SynthSeg-derived lentiform nucleus.

The lentiform nucleus was defined as the **bilateral putamen and pallidum** using SynthSeg labels:

- 12: left putamen,
- 13: left pallidum,
- 51: right putamen,
- 52: right pallidum.

A single contiguous lesion was defined as one 3D connected component using **26-connectivity**.

For the manuscript dataset, the independently reproduced selection flow was:

```text
Full cohort:                    603
Single contiguous lesion:      279
No lentiform lesion:           114
Final analysis subset:         114
```

The CST-overlap criterion was **not** used for this supplementary lentiform analysis.

The table contains the group-level values reported in Supplementary Table S8, including:

- native-space lentiform-to-lesion COM distance,
- MNI152 round-trip distance,
- EBT round-trip distance,
- MNI152 round-trip distance error,
- EBT round-trip distance error,
- between-template error difference,
- paired t-test result.

Round-trip distance error was defined as:

```text
|dist_round-trip − dist_native|
```

The subject-selection procedure can be independently reconstructed using:

```text
analysis/select_lentiform_subset_from_full_cohort.py
```

and the subject-level distances can be regenerated from a local ATLAS working directory using:

```text
analysis/compute_lentiform_lesion_dist_native_roundtrip.py
```

---

## 3.6 `data/subjects_included.xlsx`

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

Reconstructs the subject subset used for the supplementary lentiform-to-lesion COM analysis directly from the full study cohort.

The script applies the actual selection criteria used for this analysis:

1. single contiguous native-space lesion,
2. no lesion voxels within the bilateral SynthSeg-derived lentiform nucleus.

The single-lesion criterion uses 3D connected-component labeling with 26-connectivity.

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

**Main outputs**

```text
lentiform_selection_all_subjects.tsv
lentiform_final_subset.tsv
```

The first output provides an auditable subject-level record of inclusion and exclusion across the full cohort.

---

## 4.6 `analysis/compute_lentiform_lesion_dist_native_roundtrip.py`

**Purpose**

Computes the 3D Euclidean distance between the COM of the native-space lentiform nucleus and the lesion COM.

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

The manuscript Supplementary Table S8 was obtained from the single-lesion cases without lesion involvement of the lentiform nucleus.

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
S9_ventricle_fraction_dCC.tsv
```

The exact output filenames may vary slightly with the release version of the script.

Supplementary Table S8 is provided separately as:

```text
data/lentiform_lesion_dist_noLentLesion_table.tsv
```

because it derives from the dedicated lentiform-analysis pipeline described above.

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
- mask-definition sensitivity analysis.

Supplementary Table S8 is supplied directly as a released summary table and can be independently regenerated from the dedicated lentiform scripts if the required source images and intermediate files are available locally.

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

The resulting ventricular fraction can be incorporated into the released per-subject metrics table.

---

## 6.5 Lentiform-analysis subject selection

```bash
python analysis/select_lentiform_subset_from_full_cohort.py \
  --root /path/to/ATLAS_workdir
```

For the manuscript dataset, the expected selection flow is:

```text
Full cohort                         : 603
Single contiguous lesion            : 279
+ no lesion within lentiform nucleus: 114
Final eligible subset               : 114
```

A successful reproduction should end with:

```text
[OK] Selection reproduced: 603 -> 279 -> 114.
```

---

## 6.6 Lentiform-to-lesion round-trip distance

```bash
python analysis/compute_lentiform_lesion_dist_native_roundtrip.py \
  --root /path/to/ATLAS_workdir
```

The script computes subject-level native, MNI152 round-trip, and EBT round-trip lentiform-to-lesion COM distances.

For Supplementary Table S8, these measurements are restricted to the 114 subjects selected by the lentiform-analysis criteria described above.

---

# 7. Interpretation of the lentiform supplementary analysis

The supplementary lentiform analysis was designed as an additional check of relative lesion position after round-trip transformation.

It does **not** test CST overlap and does not use CST overlap as an inclusion criterion.

The purpose of the analysis is to assess whether switching standard spaces produces gross changes in the spatial relationship between a lesion and a stable deep gray-matter landmark.

The lentiform nucleus was therefore used as an anatomical anchor, and the native-space lentiform-to-lesion COM distance was compared with the corresponding distance after MNI152 and EBT round-trip transformation.

The final analysis subset contained 114 cases.

---

# 8. Data redistribution and source-data scope

This repository does **not** redistribute ATLAS R2.0 subject-level MRI images or lesion masks.

Only the following are provided:

- derived numerical metrics,
- released summary tables,
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

https://doi.org/10.5281/zenodo.19472369

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