[![DOI (concept)](https://zenodo.org/badge/DOI/[CONCEPT-DOI].svg)](https://doi.org/[CONCEPT-DOI])

# 1. Overview

This repository contains released derived data and reproducibility scripts related to our study,

*Standard-space selection in stroke magnetic resonance imaging: an elderly brain template reduces deformation bias and preserves lesion geometry*

This repository is intentionally separate from the main EBT repository, **`elderly-brain-template`**, which distributes the elderly brain template itself. The present repository focuses on the **stroke-specific comparison between EBT and MNI152** and its reproducibility package.

Specifically, it includes:

- released per-subject whole-brain registration metrics for **n = 603** included 3T ATLAS R2.0 cases,
- released lesion round-trip geometry metrics in native space,
- released lesion-wise log-Jacobian summaries and lesion volumes,
- a fairness sensitivity table based on subject-derived masks warped to each template space,
- the included subject list,
- scripts to regenerate the main manuscript tables and figures from the released TSV files.

---

# 2. Repository structure

```text
ebt-mni152-stroke-standardspace/
├─ analysis/
│  ├─ build_release_metrics_tsv.py
│  ├─ compute_lesion_logjac_volume_T1_MNI_EBT.py
│  ├─ recompute_lesion_roundtrip_dice_assd.py
│  └─ make_all_from_release.R
├─ data/
│  ├─ bg_CC_gradCC_subjectmask_warp_table.tsv
│  ├─ lesion_logjac_volume_T1_MNI_EBT.tsv
│  ├─ lesion_roundtrip_dice_assd_native.tsv
│  ├─ release_per_subject_metrics.tsv
│  ├─ S5_lesionVol_regression.tsv
│  └─ subjects_included.xlsx
└─ outputs/   # generated locally; typically not version-controlled
```

---

# 3. File contents

## 3.1 `data/release_per_subject_metrics.tsv`

**Description**  
Released per-subject whole-brain metrics for the **603 included 3T ATLAS R2.0 cases**.

**Main contents**
- background CC
- Gradient-CC
- whole-brain log-Jacobian summaries
- lesion volume
- ventricle fraction
- scanner manufacturer and image-resolution metadata

## 3.2 `data/lesion_roundtrip_dice_assd_native.tsv`

**Description**  
Released per-subject lesion round-trip geometry metrics in native space.

**Main contents**
- Dice for the MNI152 and EBT conditions
- ASSD (mm) for the MNI152 and EBT conditions

## 3.3 `data/lesion_logjac_volume_T1_MNI_EBT.tsv`

**Description**  
Released lesion-wise deformation summaries and lesion volumes.

**Main contents**
- lesion volume in native and template spaces
- lesion-wise mean/median/p5/p95 log-Jacobian for MNI152 and EBT

## 3.4 `data/bg_CC_gradCC_subjectmask_warp_table.tsv`

**Description**  
Supplementary fairness sensitivity table based on subject-derived masks (HD-BET in native space) warped to each template space.

## 3.5 `data/S5_lesionVol_regression.tsv`

**Description**  
Derived supplementary regression table summarizing lesion-volume effects on whole-brain metric differences.

**Remarks**  
This table can also be regenerated from `release_per_subject_metrics.tsv`.

## 3.6 `data/subjects_included.xlsx`

**Description**  
List of the included ATLAS R2.0 3T subjects used in the final analysis cohort.

---

# 4. Analysis scripts

## 4.1 `analysis/build_release_metrics_tsv.py`

**Purpose**  
Generates the released per-subject whole-brain metrics table from a local ATLAS-based working directory.

**Main output**
- `release_per_subject_metrics.tsv`

## 4.2 `analysis/recompute_lesion_roundtrip_dice_assd.py`

**Purpose**  
Recomputes lesion round-trip Dice and ASSD in native space using template-specific forward and inverse transforms.

**Main output**
- `lesion_roundtrip_dice_assd_native.tsv`

## 4.3 `analysis/compute_lesion_logjac_volume_T1_MNI_EBT.py`

**Purpose**  
Computes lesion-wise log-Jacobian summaries and lesion volumes for MNI152 and EBT conditions.

**Main output**
- `lesion_logjac_volume_T1_MNI_EBT.tsv`

## 4.4 `analysis/make_all_from_release.R`

**Purpose**  
Reproduces the main manuscript tables and figures, as well as key supplementary tables, from the released TSV files in `data/`.

**Main outputs**
- `Table2.tsv`
- `Table3.tsv`
- `Fig2A_CC_vs_lesionVol.pdf`
- `Fig2B_wholebrain_mean_logJ.pdf`
- `Fig3A_dDice_margins.pdf`
- `Fig3B_dASSD_margins.pdf`
- `Fig3C_lesion_mean_logJ.pdf`
- `Fig4_dCC_vs_ventricle_fraction.pdf`
- supplementary tables (`S4`–`S9`, when applicable)

---

# 5. Reproducing the main manuscript tables and figures

From the repository root, run:

```bash
Rscript analysis/make_all_from_release.R --root . --outdir outputs
```

This will generate:

- `Table2.tsv`
- `Table3.tsv`
- `Fig2A_CC_vs_lesionVol.pdf`
- `Fig2B_wholebrain_mean_logJ.pdf`
- `Fig3A_dDice_margins.pdf`
- `Fig3B_dASSD_margins.pdf`
- `Fig3C_lesion_mean_logJ.pdf`
- `Fig4_dCC_vs_ventricle_fraction.pdf`
- supplementary tables (`S4`–`S9`, when applicable)

---

# 6. Regenerating the released TSV files from a local ATLAS-derived working directory

These scripts are provided for transparency and require access to a local ATLAS-based working directory and the relevant template files.

## 6.1 Whole-brain metrics

```bash
python analysis/build_release_metrics_tsv.py --root /path/to/ATLAS_workdir
```

## 6.2 Lesion round-trip Dice/ASSD

```bash
python analysis/recompute_lesion_roundtrip_dice_assd.py --root /path/to/ATLAS_workdir
```

## 6.3 Lesion-wise log-Jacobian summaries

```bash
python analysis/compute_lesion_logjac_volume_T1_MNI_EBT.py \
  --root /path/to/ATLAS_workdir \
  --mni-t1 /path/to/MNI152_T1_1mm_brain.nii.gz \
  --ebt-t1 /path/to/Elderly_brain_T1_template_1mm.nii.gz
```

---

# 7. Data redistribution and source-data scope

This repository does **not** redistribute ATLAS R2.0 subject-level images or lesion masks.  
Only derived metrics, supplementary tables, and supporting subject lists are provided here.

Released files in `data/` are derived from the **ATLAS R2.0 Training Dataset** and are provided only for reproducibility of the present study. The included subject list is provided to document the final analysis cohort and does not replace access to the original ATLAS R2.0 dataset.

Users who wish to regenerate the released TSV files from source data must obtain ATLAS R2.0 independently and comply with its original terms of use.

ATLAS R2.0 should be cited as:

Liew SL, Lo BP, Donnelly MR, Zavaliangos-Petropulu A, Jeong JN, Barisano G, et al.  
*A large, curated, open-source stroke neuroimaging dataset to improve lesion segmentation algorithms.*  
**Scientific Data.** 2022;9:320.  
https://doi.org/10.1038/s41597-022-01401-7

---

# 8. Software environment

The released scripts were developed using:

- Python 3.11
- R 4.4
- NumPy / SciPy / pandas / NiBabel
- ANTs / ANTsPy
- FSL

Please refer to the manuscript Supplementary Table S1 for a fuller software list and versions.

---

# 9. License

## 9.1 Code

All code in `analysis/` is released under the **MIT License**.

### MIT License

Copyright (c) 2026 Kazumichi Ota

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the “Software”), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## 9.2 Released derived data

Released derived data in `data/` are provided under **CC BY 4.0**, unless otherwise noted.

You are free to:
- share
- adapt

under the terms of the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license, provided that appropriate credit is given.

## 9.3 Third-party source data

ATLAS R2.0 subject-level images and lesion masks are **not** redistributed in this repository.  
Use of source data remains subject to the original ATLAS R2.0 terms.

---

# 10. Citation

If you use these released data or scripts, please cite the associated article and the repository DOI.

**Article (provisional; update when finalized)**  
Ota K, Nakazato Y, Oyama G, *et al.*  
*Standard-space selection in stroke magnetic resonance imaging: an elderly brain template reduces deformation bias and preserves lesion geometry.*

**Repository DOI (this version)**  
`[VERSION-DOI-HERE]`

**Repository DOI (concept DOI; always latest)**  
`[CONCEPT-DOI-HERE]`

**Related template repository**  
`elderly-brain-template`

---

# 11. Related repository

- **elderly-brain-template** — main repository for the elderly brain template (EBT) itself

---

# 12. Contact

- **Corresponding Author**: Kazumichi Ota  
- **Email**: kota24@saitama-med.ac.jp

If you have any questions, feedback, or encounter issues with these files, please feel free to contact us.