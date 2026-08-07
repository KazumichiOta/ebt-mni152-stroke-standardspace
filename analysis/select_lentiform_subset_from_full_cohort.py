#!/usr/bin/env python3
"""
select_lentiform_subset_from_full_cohort.py

Purpose
-------
Reconstruct, directly from the full ATLAS_2_simple cohort, the subject
subset used for the supplementary lentiform-to-lesion round-trip
center-of-mass analysis.

Selection criteria
------------------
1. The native-space lesion consists of a single contiguous 3D component.
2. No native-space lesion voxels are present within the SynthSeg-derived
   lentiform nucleus (bilateral putamen + pallidum).

Definitions reproduced from the original analysis scripts
---------------------------------------------------------
- Single lesion:
    26-connectivity in 3D.
    Lesions with <10 voxels are classified as TOO_SMALL.
- Lentiform nucleus:
    SynthSeg labels 12, 13, 51, 52
    = bilateral putamen + pallidum.

Expected input files
--------------------
ROOT/sub-*/anat/T1lesion_mask.nii.gz
ROOT/sub-*/anat/synthseg_native/labels.nii.gz

Outputs
-------
ROOT/_qc/lentiform_selection_all_subjects.tsv
ROOT/_qc/lentiform_final_subset.tsv

Expected selection flow for the manuscript dataset
---------------------------------------------------
Full cohort:                  n = 603
Single contiguous lesion:    n = 279
No lentiform lesion:         n = 114

This script performs subject selection only. It does not calculate the
lentiform-to-lesion round-trip distances used in Supplementary Table S8.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import label


# ============================================================
# Constants reproduced from original scripts
# ============================================================

MIN_LESION_VOX = 10

# SynthSeg / FreeSurfer-style labels
# Putamen: 12 (L), 51 (R)
# Pallidum: 13 (L), 52 (R)
LENTIFORM_LABELS = [12, 13, 51, 52]

EXPECTED_TOTAL_N = 603
EXPECTED_SINGLE_N = 279
EXPECTED_FINAL_N = 114


def load_binary_mask(path: Path) -> tuple[nib.Nifti1Image, np.ndarray]:
    """Load a NIfTI image and return the image and a boolean mask."""
    img = nib.load(str(path))
    data = img.get_fdata() > 0
    return img, data


def count_components(mask: np.ndarray) -> int:
    """
    Count 3D connected components using 26-connectivity,
    reproducing the original lesion-component analysis.
    """
    structure = np.ones((3, 3, 3), dtype=np.int8)
    _, n_comp = label(mask, structure=structure)
    return int(n_comp)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct the single-lesion, lesion-free-lentiform subset "
            "used for the supplementary lentiform-to-lesion COM analysis."
        )
    )

    parser.add_argument(
        "--root",
        type=Path,
        default=Path("/Users/kazumichiota/Desktop/ATLAS_2_simple"),
        help="Root directory containing sub-* subject folders.",
    )

    parser.add_argument(
        "--expected-total-n",
        type=int,
        default=EXPECTED_TOTAL_N,
        help="Expected full cohort size. Use 0 to disable check.",
    )

    parser.add_argument(
        "--expected-single-n",
        type=int,
        default=EXPECTED_SINGLE_N,
        help="Expected number of single-lesion cases. Use 0 to disable check.",
    )

    parser.add_argument(
        "--expected-final-n",
        type=int,
        default=EXPECTED_FINAL_N,
        help="Expected final subset size. Use 0 to disable check.",
    )

    args = parser.parse_args()

    root = args.root.expanduser().resolve()

    out_dir = root / "_qc"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_all = out_dir / "lentiform_selection_all_subjects.tsv"
    out_final = out_dir / "lentiform_final_subset.tsv"

    # ========================================================
    # Input checks
    # ========================================================

    if not root.exists():
        raise FileNotFoundError(f"ROOT not found: {root}")

    subject_dirs = sorted(
        p for p in root.glob("sub-*")
        if p.is_dir()
    )

    n_total = len(subject_dirs)

    print("============================================================")
    print("Lentiform subset reconstruction")
    print("============================================================")
    print(f"[INFO] ROOT: {root}")
    print(f"[INFO] Full subject count: {n_total}")
    print()

    # ========================================================
    # Process subjects
    # ========================================================

    rows = []

    for i, sub_dir in enumerate(subject_dirs, start=1):

        sid = sub_dir.name
        anat = sub_dir / "anat"

        lesion_path = anat / "T1lesion_mask.nii.gz"
        synthseg_path = anat / "synthseg_native" / "labels.nii.gz"

        row = {
            "subject": sid,
            "native_lesion_vox": np.nan,
            "n_components": np.nan,
            "single_contiguous_lesion": False,
            "lentiform_lesion_vox": np.nan,
            "no_lentiform_lesion": False,
            "final_include": False,
            "exclusion_reason": "",
        }

        print(f"[{i:3d}/{n_total:3d}] {sid}")

        # ----------------------------------------------------
        # Required files
        # ----------------------------------------------------

        if not lesion_path.exists():
            row["exclusion_reason"] = "missing_native_lesion"
            rows.append(row)
            print("    EXCLUDE: native lesion mask missing")
            continue

        if not synthseg_path.exists():
            row["exclusion_reason"] = "missing_synthseg_native"
            rows.append(row)
            print("    EXCLUDE: native SynthSeg labels missing")
            continue

        # ====================================================
        # Criterion 1: single contiguous native lesion
        # ====================================================

        lesion_img, lesion_mask = load_binary_mask(lesion_path)

        n_vox = int(lesion_mask.sum())
        row["native_lesion_vox"] = n_vox

        if n_vox < MIN_LESION_VOX:
            row["n_components"] = 0
            row["exclusion_reason"] = "lesion_too_small"
            rows.append(row)

            print(
                f"    EXCLUDE: lesion < {MIN_LESION_VOX} voxels "
                f"(n={n_vox})"
            )
            continue

        n_comp = count_components(lesion_mask)
        row["n_components"] = n_comp

        if n_comp != 1:
            row["exclusion_reason"] = "multiple_lesion_components"
            rows.append(row)
            print(f"    EXCLUDE: n_components={n_comp}")
            continue

        row["single_contiguous_lesion"] = True

        # ====================================================
        # Criterion 2: no lesion within lentiform nucleus
        # ====================================================

        synthseg_img = nib.load(str(synthseg_path))
        synthseg_data = synthseg_img.get_fdata()

        if lesion_img.shape != synthseg_img.shape:
            row["exclusion_reason"] = "lesion_SynthSeg_shape_mismatch"
            rows.append(row)

            print(
                "    EXCLUDE: lesion/SynthSeg shape mismatch "
                f"{lesion_img.shape} vs {synthseg_img.shape}"
            )
            continue

        lentiform_mask = np.isin(
            synthseg_data,
            LENTIFORM_LABELS,
        )

        if not np.any(lentiform_mask):
            row["exclusion_reason"] = "lentiform_mask_empty"
            rows.append(row)
            print("    EXCLUDE: lentiform mask empty")
            continue

        lentiform_lesion_vox = int(
            np.logical_and(lesion_mask, lentiform_mask).sum()
        )

        row["lentiform_lesion_vox"] = lentiform_lesion_vox

        if lentiform_lesion_vox > 0:
            row["exclusion_reason"] = "lesion_within_lentiform"
            rows.append(row)

            print(
                "    EXCLUDE: lesion overlaps lentiform nucleus "
                f"({lentiform_lesion_vox} voxels)"
            )
            continue

        row["no_lentiform_lesion"] = True

        # ====================================================
        # Final inclusion
        # ====================================================

        row["final_include"] = True
        row["exclusion_reason"] = ""

        rows.append(row)

        print("    INCLUDE | single lesion | lentiform lesion=0")

    # ========================================================
    # Convert to DataFrame
    # ========================================================

    df = pd.DataFrame(rows)

    df = df.sort_values("subject").reset_index(drop=True)

    final_df = df.loc[
        df["final_include"] == True
    ].copy()

    # ========================================================
    # Selection-flow counts
    # ========================================================

    n_full = len(df)

    n_single = int(
        df["single_contiguous_lesion"]
        .fillna(False)
        .sum()
    )

    n_final = len(final_df)

    print()
    print("============================================================")
    print("Selection flow")
    print("============================================================")
    print(f"Full cohort                         : {n_full}")
    print(f"Single contiguous lesion            : {n_single}")
    print(f"+ no lesion within lentiform nucleus: {n_final}")
    print(f"Final eligible subset               : {n_final}")
    print("============================================================")
    print()

    # ========================================================
    # Write outputs
    # ========================================================

    df.to_csv(
        out_all,
        sep="\t",
        index=False,
    )

    final_columns = [
        "subject",
        "native_lesion_vox",
        "n_components",
        "single_contiguous_lesion",
        "lentiform_lesion_vox",
        "no_lentiform_lesion",
    ]

    final_df[final_columns].to_csv(
        out_final,
        sep="\t",
        index=False,
    )

    print(f"[INFO] Full selection audit: {out_all}")
    print(f"[INFO] Final subset:        {out_final}")

    # ========================================================
    # Reproduction checks
    # ========================================================

    failed = False

    if args.expected_total_n > 0 and n_full != args.expected_total_n:
        print(
            f"[ERROR] Expected full cohort n={args.expected_total_n}, "
            f"but observed n={n_full}.",
            file=sys.stderr,
        )
        failed = True

    if args.expected_single_n > 0 and n_single != args.expected_single_n:
        print(
            f"[ERROR] Expected SINGLE n={args.expected_single_n}, "
            f"but observed n={n_single}.",
            file=sys.stderr,
        )
        failed = True

    if args.expected_final_n > 0 and n_final != args.expected_final_n:
        print(
            f"[ERROR] Expected final n={args.expected_final_n}, "
            f"but observed n={n_final}.",
            file=sys.stderr,
        )
        failed = True

    if failed:
        print()
        print(
            "[ERROR] The reconstructed selection does not reproduce "
            "the manuscript dataset.",
            file=sys.stderr,
        )
        print(
            "[ERROR] Check whether the input files correspond to the "
            "same analysis version used for the manuscript.",
            file=sys.stderr,
        )
        sys.exit(1)

    print()
    print(
        f"[OK] Selection reproduced: "
        f"{n_full} -> {n_single} -> {n_final}."
    )


if __name__ == "__main__":
    main()