#!/usr/bin/env python3
"""
compute_ventricle_brain_fraction.py

Compute ventricular volume and ventricular fraction from native-space
SynthSeg segmentations for an ATLAS-derived working directory.

For each subject:

    ICV = total volume of voxels with a nonzero SynthSeg label

    ventricular volume = volume of SynthSeg labels
        4, 5, 14, 15, 43, 44

    ventricular fraction = ventricular volume / ICV

An auxiliary non-CSF-labeled volume is also retained for compatibility with
the original analysis workflow. It is defined as ICV excluding SynthSeg
label 24.

Expected input structure
------------------------

ROOT/
├── sub-XXX/
│   └── anat/
│       └── synthseg_native_wholehead_robust/
│           └── labels.nii.gz
├── sub-YYY/
│   └── anat/
│       └── synthseg_native_wholehead_robust/
│           └── labels.nii.gz
└── ...

Default output
--------------

ROOT/_qc/ventricle_brain_fraction.tsv

Dependencies
------------

numpy
nibabel

Example
-------

python analysis/compute_ventricle_brain_fraction.py \
    --root /path/to/ATLAS_workdir

Optional arguments
------------------

--output /path/to/output.tsv
--jobs 10
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import nibabel as nib
import numpy as np


# -------------------------------------------------------------------------
# SynthSeg / FreeSurfer-style label definitions
# -------------------------------------------------------------------------

# Lateral ventricles:
#   4  = left lateral ventricle
#   43 = right lateral ventricle
#
# Inferior lateral ventricles:
#   5  = left inferior lateral ventricle
#   44 = right inferior lateral ventricle
#
# Midline ventricles:
#   14 = third ventricle
#   15 = fourth ventricle
VENTRICLE_LABELS = (4, 5, 14, 15, 43, 44)

# SynthSeg / aseg CSF label.
#
# This is used only for the auxiliary non-CSF-labeled volume retained from
# the original workflow. It is not used in the ventricular-fraction
# calculation.
CSF_LABELS = (24,)

SYNTHSEG_NATIVE = os.path.join(
    "synthseg_native_wholehead_robust",
    "labels.nii.gz",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute ventricular volume and ventricular fraction from "
            "native-space SynthSeg labels."
        )
    )

    parser.add_argument(
        "--root",
        required=True,
        help=(
            "Root ATLAS-derived working directory containing "
            "sub-*/anat/synthseg_native_wholehead_robust/labels.nii.gz."
        ),
    )

    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output TSV path. "
            "Default: <root>/_qc/ventricle_brain_fraction.tsv"
        ),
    )

    parser.add_argument(
        "--jobs",
        type=int,
        default=10,
        help="Number of parallel worker processes. Default: 10.",
    )

    return parser.parse_args()


def load_labels(path: str) -> tuple[np.ndarray, float]:
    """
    Load a SynthSeg label image.

    Returns
    -------
    labels : np.ndarray
        Integer-valued 3D label array.

    voxel_volume_mm3 : float
        Physical voxel volume in mm^3.
    """
    img = nib.load(path)

    labels = np.asarray(
        np.rint(img.get_fdata()),
        dtype=np.int32,
    )

    zooms = img.header.get_zooms()[:3]
    voxel_volume_mm3 = float(
        zooms[0] * zooms[1] * zooms[2]
    )

    return labels, voxel_volume_mm3


def process_subject(
    root: str,
    subject: str,
) -> tuple:
    """
    Process one subject.

    Returns
    -------
    ("OK", subject, row_dict)

    or

    ("SKIP", subject, reason)
    """
    anat_dir = os.path.join(
        root,
        subject,
        "anat",
    )

    label_path = os.path.join(
        anat_dir,
        SYNTHSEG_NATIVE,
    )

    if not os.path.exists(label_path):
        return (
            "SKIP",
            subject,
            f"missing {label_path}",
        )

    try:
        labels, voxel_volume_mm3 = load_labels(
            label_path
        )
    except Exception as exc:
        return (
            "SKIP",
            subject,
            f"failed to load labels: {exc}",
        )

    if labels.ndim != 3:
        return (
            "SKIP",
            subject,
            f"expected 3D label image, got shape {labels.shape}",
        )

    # ------------------------------------------------------------------
    # Intracranial volume
    #
    # Definition used in the manuscript analysis:
    # all voxels assigned a nonzero SynthSeg label.
    # ------------------------------------------------------------------
    icv_mask = labels > 0

    icv_voxels = int(
        np.count_nonzero(icv_mask)
    )

    if icv_voxels == 0:
        return (
            "SKIP",
            subject,
            "ICV is zero (no nonzero SynthSeg labels)",
        )

    # ------------------------------------------------------------------
    # Ventricular volume
    # ------------------------------------------------------------------
    ventricle_mask = np.isin(
        labels,
        VENTRICLE_LABELS,
    )

    ventricle_voxels = int(
        np.count_nonzero(ventricle_mask)
    )

    # ------------------------------------------------------------------
    # Auxiliary measure retained from the original workflow.
    #
    # This is ICV excluding label 24 only. It is not required for the
    # ventricular-fraction analysis and should not be interpreted as a
    # strict tissue-segmented parenchymal volume.
    # ------------------------------------------------------------------
    csf_mask = np.isin(
        labels,
        CSF_LABELS,
    )

    non_csf_mask = (
        icv_mask
        & (~csf_mask)
    )

    non_csf_voxels = int(
        np.count_nonzero(non_csf_mask)
    )

    # ------------------------------------------------------------------
    # Physical volumes
    # ------------------------------------------------------------------
    icv_mm3 = (
        icv_voxels
        * voxel_volume_mm3
    )

    ventricle_mm3 = (
        ventricle_voxels
        * voxel_volume_mm3
    )

    brain_mm3 = (
        non_csf_voxels
        * voxel_volume_mm3
    )

    ventricle_fraction = (
        ventricle_mm3
        / icv_mm3
    )

    brain_fraction = (
        brain_mm3
        / icv_mm3
    )

    row = {
        "subject": subject,
        "ICV_mm3": icv_mm3,
        "ventricle_mm3": ventricle_mm3,
        "brain_mm3": brain_mm3,
        "ventricle_fraction": ventricle_fraction,
        "brain_fraction": brain_fraction,
    }

    return (
        "OK",
        subject,
        row,
    )


def main() -> None:
    args = parse_args()

    root = os.path.abspath(
        os.path.expanduser(args.root)
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
            "ventricle_brain_fraction.tsv",
        )
    else:
        output_tsv = os.path.abspath(
            os.path.expanduser(args.output)
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
    # Discover subject directories
    # ------------------------------------------------------------------
    subjects = sorted(
        os.path.basename(path)
        for path in glob.glob(
            os.path.join(
                root,
                "sub-*",
            )
        )
        if os.path.isdir(path)
    )

    n_total = len(subjects)

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

    results_ok = []

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
            for subject in subjects
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

            percentage = (
                100.0
                * index
                / n_total
            )

            try:
                result = future.result()
            except Exception as exc:
                n_skip += 1

                print(
                    f"[{index}/{n_total} "
                    f"{percentage:5.1f}%] "
                    f"{subject}: "
                    f"EXCEPTION ({exc})"
                )
                continue

            status = result[0]

            if status == "OK":
                _, subject, row = result

                results_ok.append(row)
                n_ok += 1

                print(
                    f"[{index}/{n_total} "
                    f"{percentage:5.1f}%] "
                    f"{subject}: OK "
                    f"(ventricle_fraction="
                    f"{row['ventricle_fraction']:.6f})"
                )

            else:
                _, subject, reason = result

                n_skip += 1

                print(
                    f"[{index}/{n_total} "
                    f"{percentage:5.1f}%] "
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
    # Write TSV
    # ------------------------------------------------------------------
    results_ok.sort(
        key=lambda row: row["subject"]
    )

    fieldnames = [
        "subject",
        "ICV_mm3",
        "ventricle_mm3",
        "brain_mm3",
        "ventricle_fraction",
        "brain_fraction",
    ]

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
            fieldnames
        )

        for row in results_ok:
            writer.writerow(
                [
                    row["subject"],
                    f"{row['ICV_mm3']:.2f}",
                    f"{row['ventricle_mm3']:.2f}",
                    f"{row['brain_mm3']:.2f}",
                    f"{row['ventricle_fraction']:.6f}",
                    f"{row['brain_fraction']:.6f}",
                ]
            )

    print(
        f"[INFO] output written: "
        f"{output_tsv}"
    )


if __name__ == "__main__":
    main()