#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
recompute_lesion_roundtrip_dice_assd.py

Recompute lesion round-trip Dice and ASSD for all subjects:
  - native lesion mask: anat/T1lesion_mask.nii.gz
  - MNI roundtrip mask: anat/T1lesion_MNI_roundtrip.nii.gz
  - EBT roundtrip mask: anat/T1lesion_EBT_roundtrip.nii.gz

ASSD computed in native space using surface distances (mm),
with voxel spacing taken from native lesion NIfTI header.
Outputs:
  ROOT/_qc/lesion_roundtrip_dice_assd_native.tsv
"""

import argparse, os, glob
import numpy as np
import nibabel as nib
from scipy.ndimage import binary_erosion, distance_transform_edt

def load_bool(path):
    img = nib.load(path)
    data = img.get_fdata()
    return (data > 0.5), img

def dice(a, b):
    a = a.astype(bool); b = b.astype(bool)
    inter = np.logical_and(a,b).sum()
    denom = a.sum() + b.sum()
    if denom == 0:
        return np.nan
    return float(2.0 * inter / denom)

def surface(mask):
    # boundary voxels: mask XOR eroded(mask)
    if mask.sum() == 0:
        return mask
    er = binary_erosion(mask, structure=np.ones((3,3,3), dtype=bool), iterations=1)
    return np.logical_and(mask, np.logical_not(er))

def assd(a, b, spacing):
    # average symmetric surface distance (mm)
    sa = surface(a); sb = surface(b)
    if sa.sum()==0 or sb.sum()==0:
        return np.nan
    # distance to surface B
    dt_b = distance_transform_edt(~sb, sampling=spacing)
    dt_a = distance_transform_edt(~sa, sampling=spacing)
    d_ab = dt_b[sa].mean()
    d_ba = dt_a[sb].mean()
    return float((d_ab + d_ba)/2.0)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/Volumes/Extreme Pro/ATLAS2.0/ATLAS_2_simple")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = args.root
    out_path = args.out or os.path.join(root, "_qc", "lesion_roundtrip_dice_assd_native.tsv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    subs = sorted([os.path.basename(p) for p in glob.glob(os.path.join(root, "sub-*")) if os.path.isdir(p)])
    if not subs:
        raise RuntimeError("No sub-* under root")

    rows = []
    skip = 0
    for sid in subs:
        anat = os.path.join(root, sid, "anat")
        nat_p = os.path.join(anat, "T1lesion_mask.nii.gz")
        mni_p = os.path.join(anat, "T1lesion_MNI_roundtrip.nii.gz")
        ebt_p = os.path.join(anat, "T1lesion_EBT_roundtrip.nii.gz")
        if not (os.path.exists(nat_p) and os.path.exists(mni_p) and os.path.exists(ebt_p)):
            skip += 1
            continue

        nat, nat_img = load_bool(nat_p)
        mni, _ = load_bool(mni_p)
        ebt, _ = load_bool(ebt_p)

        if nat.shape != mni.shape or nat.shape != ebt.shape:
            skip += 1
            continue

        spacing = nat_img.header.get_zooms()[:3]
        d_mni = dice(nat, mni)
        d_ebt = dice(nat, ebt)
        a_mni = assd(nat, mni, spacing)
        a_ebt = assd(nat, ebt, spacing)

        rows.append([sid, d_mni, d_ebt, a_mni, a_ebt])

    import pandas as pd
    df = pd.DataFrame(rows, columns=["subject","dice_MNI","dice_EBT","assd_MNI_mm","assd_EBT_mm"])
    df.to_csv(out_path, sep="\t", index=False)
    print("[INFO] wrote:", out_path)
    print("[INFO] n:", df.shape[0], "skip:", skip)

if __name__ == "__main__":
    main()
