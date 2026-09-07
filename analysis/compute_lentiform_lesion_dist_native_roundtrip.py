#!/usr/bin/env python3
"""
compute_lentiform_lesion_dist_native_roundtrip.py

Compute native-space lentiform-to-lesion center-of-mass (COM) distances
for the original lesion mask and for lesion masks after round-trip
transformation through MNI152 or the elderly brain template (EBT).

The lentiform nucleus is defined as the bilateral putamen and pallidum
using SynthSeg / FreeSurfer-style labels:

    12 = left putamen
    13 = left pallidum
    51 = right putamen
    52 = right pallidum

For each subject, the script calculates:

    dist_native_mm
        Distance between the lentiform COM and the native lesion COM.

    dist_MNIround_mm
        Distance between the lentiform COM and the lesion COM after
        native -> MNI152 -> native round-trip transformation.

    dist_EBTround_mm
        Distance between the lentiform COM and the lesion COM after
        native -> EBT -> native round-trip transformation.

    err_MNI_mm
        |dist_MNIround_mm - dist_native_mm|

    err_EBT_mm
        |dist_EBTround_mm - dist_native_mm|

Expected input structure
------------------------

ROOT/
├── sub-XXX/
│   └── anat/
│       ├── T1lesion_mask.nii.gz
│       ├── T1lesion_MNI_roundtrip.nii.gz
│       ├── T1lesion_EBT_roundtrip.nii.gz
│       └── synthseg_native_wholehead_robust/
│           └── labels.nii.gz
└── ...

The SynthSeg label maps used in the final analysis were generated in
robust mode from N4-corrected native whole-head T1-weighted images before
skull stripping.

Default output
--------------

ROOT/_qc/lentiform_lesion_dist_native_roundtrip.tsv

Important coordinate convention
-------------------------------

NumPy indices obtained from a NiBabel image are voxel coordinates in
(i, j, k) order. These coordinates are passed directly to the NIfTI affine
using nibabel.affines.apply_affine().

No axis reversal is performed.

Dependencies
------------

numpy
nibabel

Example
-------

python analysis/compute_lentiform_lesion_dist_native_roundtrip.py \
    --root /path/to/ATLAS_workdir

Optional arguments
------------------

--output /path/to/output.tsv
--jobs 5
"""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import nibabel as nib
import numpy as np


# -------------------------------------------------------------------------
# Input filenames
# -------------------------------------------------------------------------

SYNTHSEG_LABELS = os.path.join(
    "synthseg_native_wholehead_robust",
    "labels.nii.gz",
)

LESION_NATIVE = (
    "T1lesion_mask.nii.gz"
)

LESION_MNI_ROUNDTRIP = (
    "T1lesion_MNI_roundtrip.nii.gz"
)

LESION_EBT_ROUNDTRIP = (
    "T1lesion_EBT_roundtrip.nii.gz"
)


# -------------------------------------------------------------------------
# SynthSeg / FreeSurfer-style lentiform labels
# -------------------------------------------------------------------------

LENTIFORM_LABELS = (
    12,  # left putamen
    13,  # left pallidum
    51,  # right putamen
    52,  # right pallidum
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute native-space lentiform-to-lesion COM distances "
            "for native, MNI152 round-trip, and EBT round-trip "
            "lesion masks."
        )
    )

    parser.add_argument(
        "--root",
        required=True,
        help=(
            "Root ATLAS-derived working directory containing "
            "sub-*/anat directories."
        ),
    )

    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output TSV path. Default: "
            "<root>/_qc/"
            "lentiform_lesion_dist_native_roundtrip.tsv"
        ),
    )

    parser.add_argument(
        "--jobs",
        type=int,
        default=5,
        help=(
            "Number of parallel worker processes. "
            "Default: 5."
        ),
    )

    return parser.parse_args()


def voxel_com(
    mask: np.ndarray,
) -> np.ndarray | None:
    """
    Compute the center of mass of a binary mask in voxel coordinates.

    Parameters
    ----------
    mask : np.ndarray
        Three-dimensional binary or integer mask.

    Returns
    -------
    np.ndarray or None
        Mean voxel coordinate in NiBabel-compatible (i, j, k) order.
        Returns None if the mask is empty.
    """
    indices = np.argwhere(
        mask > 0
    )

    if indices.size == 0:
        return None

    return indices.mean(
        axis=0
    )


def mask_com_world(
    image: nib.spatialimages.SpatialImage,
    mask: np.ndarray,
) -> np.ndarray | None:
    """
    Compute a mask center of mass in world coordinates (mm).

    The voxel COM returned by np.argwhere() is in array index
    coordinates (i, j, k), which correspond directly to the voxel
    coordinates expected by the NIfTI affine.

    Therefore the coordinate order must NOT be reversed.
    """
    com_ijk = voxel_com(
        mask
    )

    if com_ijk is None:
        return None

    com_world = nib.affines.apply_affine(
        image.affine,
        com_ijk,
    )

    return np.asarray(
        com_world,
        dtype=float,
    )


def lesion_com_world(
    path: str,
) -> np.ndarray | None:
    """
    Load a lesion mask and return its COM in world coordinates (mm).

    Any voxel with a value > 0 is treated as lesion.
    """
    image = nib.load(
        path
    )

    data = image.get_fdata()

    if data.ndim != 3:
        raise ValueError(
            f"Expected 3D lesion image at {path}, "
            f"got shape {data.shape}"
        )

    lesion_mask = (
        data > 0
    )

    return mask_com_world(
        image,
        lesion_mask,
    )


def euclidean_distance(
    point_a: np.ndarray,
    point_b: np.ndarray,
) -> float:
    """
    Compute 3D Euclidean distance in millimeters.
    """
    return float(
        np.linalg.norm(
            point_a - point_b
        )
    )


def process_subject(
    root: str,
    subject: str,
) -> tuple:
    """
    Process one subject.

    Returns
    -------
    ("OK", subject, dist_native, dist_MNIround, dist_EBTround,
     err_MNI, err_EBT)

    or

    ("SKIP", subject, reason)
    """
    try:
        anat_dir = os.path.join(
            root,
            subject,
            "anat",
        )

        label_path = os.path.join(
            anat_dir,
            SYNTHSEG_LABELS,
        )

        lesion_native_path = os.path.join(
            anat_dir,
            LESION_NATIVE,
        )

        lesion_mni_roundtrip_path = os.path.join(
            anat_dir,
            LESION_MNI_ROUNDTRIP,
        )

        lesion_ebt_roundtrip_path = os.path.join(
            anat_dir,
            LESION_EBT_ROUNDTRIP,
        )

        required_files = [
            label_path,
            lesion_native_path,
            lesion_mni_roundtrip_path,
            lesion_ebt_roundtrip_path,
        ]

        missing_files = [
            path
            for path in required_files
            if not os.path.exists(path)
        ]

        if missing_files:
            missing_names = [
                os.path.relpath(
                    path,
                    anat_dir,
                )
                for path in missing_files
            ]

            return (
                "SKIP",
                subject,
                "missing: "
                + ", ".join(
                    missing_names
                ),
            )

        # --------------------------------------------------------------
        # Lentiform nucleus COM
        # --------------------------------------------------------------
        label_image = nib.load(
            label_path
        )

        label_data = np.asarray(
            np.rint(
                label_image.get_fdata()
            ),
            dtype=np.int32,
        )

        if label_data.ndim != 3:
            return (
                "SKIP",
                subject,
                (
                    "SynthSeg label image is not 3D: "
                    f"{label_data.shape}"
                ),
            )

        lentiform_mask = np.isin(
            label_data,
            LENTIFORM_LABELS,
        )

        if not np.any(
            lentiform_mask
        ):
            return (
                "SKIP",
                subject,
                "lentiform mask is empty",
            )

        lentiform_com_world = (
            mask_com_world(
                label_image,
                lentiform_mask,
            )
        )

        if lentiform_com_world is None:
            return (
                "SKIP",
                subject,
                "lentiform COM could not be calculated",
            )

        # --------------------------------------------------------------
        # Lesion COMs
        # --------------------------------------------------------------
        native_com_world = lesion_com_world(
            lesion_native_path
        )

        mni_roundtrip_com_world = (
            lesion_com_world(
                lesion_mni_roundtrip_path
            )
        )

        ebt_roundtrip_com_world = (
            lesion_com_world(
                lesion_ebt_roundtrip_path
            )
        )

        if native_com_world is None:
            return (
                "SKIP",
                subject,
                "native lesion mask is empty",
            )

        if mni_roundtrip_com_world is None:
            return (
                "SKIP",
                subject,
                "MNI152 round-trip lesion mask is empty",
            )

        if ebt_roundtrip_com_world is None:
            return (
                "SKIP",
                subject,
                "EBT round-trip lesion mask is empty",
            )

        # --------------------------------------------------------------
        # Lentiform-to-lesion distances
        # --------------------------------------------------------------
        dist_native = euclidean_distance(
            lentiform_com_world,
            native_com_world,
        )

        dist_mni_roundtrip = (
            euclidean_distance(
                lentiform_com_world,
                mni_roundtrip_com_world,
            )
        )

        dist_ebt_roundtrip = (
            euclidean_distance(
                lentiform_com_world,
                ebt_roundtrip_com_world,
            )
        )

        # --------------------------------------------------------------
        # Absolute round-trip distance errors
        # --------------------------------------------------------------
        err_mni = abs(
            dist_mni_roundtrip
            - dist_native
        )

        err_ebt = abs(
            dist_ebt_roundtrip
            - dist_native
        )

        return (
            "OK",
            subject,
            dist_native,
            dist_mni_roundtrip,
            dist_ebt_roundtrip,
            err_mni,
            err_ebt,
        )

    except Exception as exc:
        return (
            "SKIP",
            subject,
            f"exception: {repr(exc)}",
        )


def main() -> None:
    args = parse_args()

    root = os.path.abspath(
        os.path.expanduser(
            args.root
        )
    )

    if not os.path.isdir(root):
        raise FileNotFoundError(
            f"Root directory not found: {root}"
        )

    if args.jobs < 1:
        raise ValueError(
            "--jobs must be >= 1"
        )

    if args.output is None:
        output_tsv = os.path.join(
            root,
            "_qc",
            "lentiform_lesion_dist_native_roundtrip.tsv",
        )
    else:
        output_tsv = os.path.abspath(
            os.path.expanduser(
                args.output
            )
        )

    output_dir = os.path.dirname(
        output_tsv
    )

    if output_dir:
        os.makedirs(
            output_dir,
            exist_ok=True,
        )

    # ------------------------------------------------------------------
    # Discover subjects
    # ------------------------------------------------------------------
    subject_ids = sorted(
        directory
        for directory in os.listdir(
            root
        )
        if (
            directory.startswith(
                "sub-"
            )
            and os.path.isdir(
                os.path.join(
                    root,
                    directory,
                )
            )
        )
    )

    n_total = len(
        subject_ids
    )

    print(
        f"[INFO] root: {root}"
    )
    print(
        f"[INFO] subjects found: {n_total}"
    )
    print(
        f"[INFO] parallel workers: {args.jobs}"
    )

    if n_total == 0:
        print(
            "[WARN] No sub-* subject directories found."
        )
        return

    results = []

    n_ok = 0
    n_skip = 0

    # ------------------------------------------------------------------
    # Parallel processing
    # ------------------------------------------------------------------
    with ProcessPoolExecutor(
        max_workers=args.jobs
    ) as executor:

        future_to_subject = {
            executor.submit(
                process_subject,
                root,
                subject,
            ): subject
            for subject in subject_ids
        }

        for index, future in enumerate(
            as_completed(
                future_to_subject
            ),
            start=1,
        ):
            subject = future_to_subject[
                future
            ]

            try:
                result = future.result()
            except Exception as exc:
                n_skip += 1

                print(
                    f"[{index}/{n_total}] "
                    f"{subject}: "
                    f"EXCEPTION ({repr(exc)})"
                )
                continue

            status = result[0]

            if status == "OK":
                (
                    _,
                    subject,
                    dist_native,
                    dist_mni_roundtrip,
                    dist_ebt_roundtrip,
                    err_mni,
                    err_ebt,
                ) = result

                results.append(
                    (
                        subject,
                        dist_native,
                        dist_mni_roundtrip,
                        dist_ebt_roundtrip,
                        err_mni,
                        err_ebt,
                    )
                )

                n_ok += 1

                print(
                    f"[{index}/{n_total}] "
                    f"{subject}: OK"
                )

            else:
                _, subject, reason = result

                n_skip += 1

                print(
                    f"[{index}/{n_total}] "
                    f"{subject}: "
                    f"SKIP ({reason})"
                )

    print(
        f"[INFO] finished: "
        f"OK={n_ok}, "
        f"SKIP={n_skip}, "
        f"TOTAL={n_total}"
    )

    if n_ok == 0:
        print(
            "[WARN] No valid subjects. "
            "Output TSV will not be written."
        )
        return

    # ------------------------------------------------------------------
    # Write subject-level TSV
    # ------------------------------------------------------------------
    results.sort(
        key=lambda row: row[0]
    )

    with open(
        output_tsv,
        "w",
        newline="",
        encoding="utf-8",
    ) as file_obj:

        writer = csv.writer(
            file_obj,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "subject",
                "dist_native_mm",
                "dist_MNIround_mm",
                "dist_EBTround_mm",
                "err_MNI_mm",
                "err_EBT_mm",
            ]
        )

        for row in results:
            writer.writerow(
                [
                    row[0],
                    f"{row[1]:.10f}",
                    f"{row[2]:.10f}",
                    f"{row[3]:.10f}",
                    f"{row[4]:.10f}",
                    f"{row[5]:.10f}",
                ]
            )

    print(
        f"[INFO] output written: "
        f"{output_tsv}"
    )


if __name__ == "__main__":
    main()