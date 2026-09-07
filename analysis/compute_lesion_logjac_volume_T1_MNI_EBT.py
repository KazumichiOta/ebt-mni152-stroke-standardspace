#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_lesion_logjac_volume_T1_MNI_EBT.py

Compute lesion volumes and intralesional log-Jacobian summary statistics
for MNI152 and EBT template spaces.

This script is intended for public release / reproducibility:
- No hard-coded user paths (all paths via CLI)
- Clear expected inputs and outputs
- Robust missing-file handling with summary
- Default behavior matches manuscript definition:
    logJ = log(detJ(affine+SyN)) excluding rigid
  implemented as:
    logJ_total = logJ_warp + log(det(Affine))   (Affine term is spatially constant)

Dependencies:
  - Python 3.8+
  - numpy
  - antspyx / ANTsPy (import ants)

Expected per-subject files under: <root>/sub-*/anat/
  Native lesion mask:
    T1lesion_mask.nii.gz
  Lesion masks already warped to template spaces (nearest-neighbor):
    T1lesion_in_MNI_syn_NN.nii.gz
    T1lesion_in_EBT_syn_NN.nii.gz
  Forward SyN warp fields (native -> template):
    T1w2MNI_syn_1Warp.nii.gz
    T1w2EBT_syn_1Warp.nii.gz
  Forward affine transforms (native -> template):
    T1w2MNI_syn_0GenericAffine.mat
    T1w2EBT_syn_0GenericAffine.mat

Outputs:
  <root>/_qc/lesion_logjac_volume_T1_MNI_EBT.tsv
    subject
    lesion_vol_T1_mm3, lesion_vol_MNI_mm3, lesion_vol_EBT_mm3
    logJac_mean/median/p5/p95_{MNI,EBT}_lesion

Usage (example):
  conda activate antspy
  python3 compute_lesion_logjac_volume_T1_MNI_EBT.py \
    --root "/path/to/ATLAS_workdir" \
    --mni-t1 "/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz" \
    --ebt-t1 "/path/to/Elderly_brain_T1_template_1mm.nii.gz"
"""

import argparse
import csv
import glob
import math
import os
import re
import sys
import time
from typing import Optional, Tuple, List

import numpy as np
import ants


# ----------------------------
# I/O helpers
# ----------------------------
def die(msg: str, code: int = 1) -> None:
    print(f"[FATAL] {msg}", file=sys.stderr)
    sys.exit(code)


def ensure_file(path: str, label: str) -> None:
    if not os.path.exists(path):
        die(f"{label} not found: {path}")


def read_mask(path: str) -> ants.ANTsImage:
    """Read mask as ANTsImage; binarization is done by threshold (>0.5)."""
    return ants.image_read(path)


def mask_volume_mm3(mask_img: ants.ANTsImage) -> float:
    """Compute volume (mm^3) from binary mask (>0.5)."""
    arr = (mask_img.numpy() > 0.5)
    vox = int(arr.sum())
    sp = mask_img.spacing
    return float(vox * sp[0] * sp[1] * sp[2])


# ----------------------------
# Affine log(det) parsing
# ----------------------------
def parse_ants_affine_3x3(mat_path: str) -> Optional[np.ndarray]:
    """
    Parse ANTs/ITK .mat file (MatrixOffsetTransformBase_double_3_3).
    Returns a 3x3 matrix (numpy) if found, otherwise None.
    """
    if not os.path.exists(mat_path):
        return None

    txt = open(mat_path, "r", encoding="utf-8", errors="replace").read()
    m = re.search(r"^Parameters:\s+(.+)$", txt, flags=re.MULTILINE)
    if not m:
        return None
    nums = [float(x) for x in m.group(1).strip().split()]
    if len(nums) < 9:
        return None
    return np.array(nums[:9], dtype=np.float64).reshape((3, 3))


def affine_logdet(mat_path: str) -> float:
    """
    Return log(det(A)) from affine matrix file.
    - Returns NaN if missing/unparseable/invalid (det<=0).
    """
    M = parse_ants_affine_3x3(mat_path)
    if M is None:
        return float("nan")
    det = float(np.linalg.det(M))
    if det <= 0:
        return float("nan")
    return float(math.log(det))


# ----------------------------
# logJac stats within lesion
# ----------------------------
def jacobian_log_image(domain_img: ants.ANTsImage, warp_path: str) -> ants.ANTsImage:
    """
    Create log-Jacobian image for a forward warp field.
    `domain_img` defines the grid (template space).
    """
    return ants.create_jacobian_determinant_image(
        domain_img,
        warp_path,
        True,   # do_log=True
        False,  # geom=False
    )


def stats_in_mask(values: np.ndarray) -> Tuple[float, float, float, float]:
    """Return mean, median, p5, p95."""
    return (
        float(np.mean(values)),
        float(np.median(values)),
        float(np.percentile(values, 5)),
        float(np.percentile(values, 95)),
    )


def lesion_logjac_stats(
    template_domain: ants.ANTsImage,
    warp_path: str,
    lesion_mask_template: ants.ANTsImage,
    affine_mat_path: Optional[str],
    include_affine: bool,
) -> Tuple[float, float, float, float]:
    """
    Compute intralesional logJac stats in template space.
    Default (include_affine=True):
      logJac_total = logJac_warp + log(det(Affine))
    """
    if not os.path.exists(warp_path):
        return (np.nan, np.nan, np.nan, np.nan)

    jac_log_img = jacobian_log_image(template_domain, warp_path)
    jac_log = jac_log_img.numpy()

    if include_affine and affine_mat_path:
        ld = affine_logdet(affine_mat_path)
        if np.isfinite(ld):
            jac_log = jac_log + ld  # constant shift

    mask = (lesion_mask_template.numpy() > 0.5)
    vals = jac_log[mask]
    if vals.size == 0:
        return (np.nan, np.nan, np.nan, np.nan)

    return stats_in_mask(vals)


# ----------------------------
# Main
# ----------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Root directory containing sub-*/anat/")
    ap.add_argument("--out", default=None, help="Output TSV (default: <root>/_qc/lesion_logjac_volume_T1_MNI_EBT.tsv)")

    ap.add_argument("--mni-t1", required=True, help="MNI template T1 (brain-only) used as Jacobian domain")
    ap.add_argument("--ebt-t1", required=True, help="EBT template T1 (brain-only) used as Jacobian domain")

    # Default: match manuscript (affine+SyN; rigid excluded)
    ap.add_argument("--no-affine", action="store_true", help="Compute warp-only logJac stats (disables affine logdet addition)")

    # File name conventions (defaults = your pipeline)
    ap.add_argument("--lesion-t1-name", default="T1lesion_mask.nii.gz")
    ap.add_argument("--lesion-mni-name", default="T1lesion_in_MNI_syn_NN.nii.gz")
    ap.add_argument("--lesion-ebt-name", default="T1lesion_in_EBT_syn_NN.nii.gz")
    ap.add_argument("--mni-warp-name", default="T1w2MNI_syn_1Warp.nii.gz")
    ap.add_argument("--ebt-warp-name", default="T1w2EBT_syn_1Warp.nii.gz")
    ap.add_argument("--mni-affine-name", default="T1w2MNI_syn_0GenericAffine.mat")
    ap.add_argument("--ebt-affine-name", default="T1w2EBT_syn_0GenericAffine.mat")

    ap.add_argument("--min-lesion-vox", type=int, default=1, help="Skip subjects with lesion voxels < this (default 1)")
    args = ap.parse_args()

    root = args.root
    out_tsv = args.out or os.path.join(root, "_qc", "lesion_logjac_volume_T1_MNI_EBT.tsv")
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)

    ensure_file(args.mni_t1, "MNI template T1")
    ensure_file(args.ebt_t1, "EBT template T1")

    include_affine = not args.no_affine

    print("[INFO] Loading templates...")
    tpl_mni = ants.image_read(args.mni_t1)
    tpl_ebt = ants.image_read(args.ebt_t1)
    print(f"  MNI domain: shape={tpl_mni.shape}, spacing={tpl_mni.spacing}")
    print(f"  EBT domain: shape={tpl_ebt.shape}, spacing={tpl_ebt.spacing}")
    print(f"[INFO] logJac mode: {'warp+affine (rigid excluded)' if include_affine else 'warp-only'}")
    print(f"[INFO] Output: {out_tsv}")
    print("")

    sub_dirs = sorted([p for p in glob.glob(os.path.join(root, "sub-*")) if os.path.isdir(p)])
    if not sub_dirs:
        die(f"No sub-* directories under: {root}")

    header = [
        "subject",
        "lesion_vol_T1_mm3",
        "lesion_vol_MNI_mm3",
        "lesion_vol_EBT_mm3",
        "logJac_mean_MNI_lesion",
        "logJac_median_MNI_lesion",
        "logJac_p5_MNI_lesion",
        "logJac_p95_MNI_lesion",
        "logJac_mean_EBT_lesion",
        "logJac_median_EBT_lesion",
        "logJac_p5_EBT_lesion",
        "logJac_p95_EBT_lesion",
    ]

    n_ok, n_skip = 0, 0
    t0 = time.time()

    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)

        for i, sub_dir in enumerate(sub_dirs, start=1):
            sid = os.path.basename(sub_dir)
            anat = os.path.join(sub_dir, "anat")
            if not os.path.isdir(anat):
                n_skip += 1
                continue

            p_t1 = os.path.join(anat, args.lesion_t1_name)
            p_mni_les = os.path.join(anat, args.lesion_mni_name)
            p_ebt_les = os.path.join(anat, args.lesion_ebt_name)
            p_mni_warp = os.path.join(anat, args.mni_warp_name)
            p_ebt_warp = os.path.join(anat, args.ebt_warp_name)
            p_mni_aff = os.path.join(anat, args.mni_affine_name)
            p_ebt_aff = os.path.join(anat, args.ebt_affine_name)

            needed = [p_t1, p_mni_les, p_ebt_les, p_mni_warp, p_ebt_warp]
            if any(not os.path.exists(p) for p in needed):
                n_skip += 1
                continue

            # Read lesion masks
            T1_les = read_mask(p_t1)
            # skip tiny lesions if requested
            if int((T1_les.numpy() > 0.5).sum()) < args.min_lesion_vox:
                n_skip += 1
                continue

            MNI_les = read_mask(p_mni_les)
            EBT_les = read_mask(p_ebt_les)

            # Volumes
            vol_t1 = mask_volume_mm3(T1_les)
            vol_mni = mask_volume_mm3(MNI_les)
            vol_ebt = mask_volume_mm3(EBT_les)

            # logJac stats in template space (within template-space lesion masks)
            mni_stats = lesion_logjac_stats(
                template_domain=tpl_mni,
                warp_path=p_mni_warp,
                lesion_mask_template=MNI_les,
                affine_mat_path=(p_mni_aff if os.path.exists(p_mni_aff) else None),
                include_affine=include_affine,
            )
            ebt_stats = lesion_logjac_stats(
                template_domain=tpl_ebt,
                warp_path=p_ebt_warp,
                lesion_mask_template=EBT_les,
                affine_mat_path=(p_ebt_aff if os.path.exists(p_ebt_aff) else None),
                include_affine=include_affine,
            )

            w.writerow([
                sid,
                vol_t1, vol_mni, vol_ebt,
                *mni_stats,
                *ebt_stats
            ])
            n_ok += 1

            if i % 25 == 0 or i == len(sub_dirs):
                print(f"[INFO] {i}/{len(sub_dirs)} processed (OK={n_ok}, SKIP={n_skip})")

    dt = time.time() - t0
    print("")
    print("[INFO] Done.")
    print("[INFO] OK:", n_ok, "SKIP:", n_skip)
    print("[INFO] Elapsed (min):", round(dt / 60.0, 2))


if __name__ == "__main__":
    main()
