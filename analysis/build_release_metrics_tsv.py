#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
build_release_metrics_tsv.py

Creates a release-ready per-subject TSV by merging:
  1) Harmonized CC/gradient-CC (brain-only same-rule masks)
     ROOT/_qc/bg_CC_gradCC_brainonly_harmonized.tsv
  2) Whole-brain logJac metrics + covariates (IGNORE old CC/SSIM/COM columns)
     ROOT/metrics_ATLAS_to_MNI_EBT_with_covariates.tsv
  3) Lesion volume (native) for Fig2A x-axis
     ROOT/_qc/lesion_fraction_qc_v2.tsv (case, Vlesion_native_mm3)
  4) Atrophy proxy
     ROOT/_qc/ventricle_brain_fraction.tsv (ventricle_fraction, brain_fraction)

Outputs:
  ROOT/_qc/release_per_subject_metrics.tsv
"""

import argparse
import os
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/Volumes/Extreme Pro/ATLAS2.0/ATLAS_2_simple")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = args.root
    out_path = args.out or os.path.join(root, "_qc", "release_per_subject_metrics.tsv")

    cc_path = os.path.join(root, "_qc", "bg_CC_gradCC_brainonly_harmonized.tsv")
    logj_path = os.path.join(root, "metrics_ATLAS_to_MNI_EBT_with_covariates.tsv")
    lesion_path = os.path.join(root, "_qc", "lesion_fraction_qc_v2.tsv")
    vf_path = os.path.join(root, "_qc", "ventricle_brain_fraction.tsv")

    for p in [cc_path, logj_path, lesion_path, vf_path]:
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    cc = pd.read_csv(cc_path, sep="\t")
    # expected cols: subject, CC_MNI_bg, CC_EBT_bg, gCC_MNI_bg, gCC_EBT_bg (and maybe others)
    need_cc = ["subject", "CC_MNI_bg", "CC_EBT_bg"]
    if not all(c in cc.columns for c in need_cc):
        raise RuntimeError(f"Missing columns in CC TSV. Have: {list(cc.columns)}")

    # logJ/covariates TSV: keep ONLY logJac + covariates (drop old CC/SSIM/COM)
    lj = pd.read_csv(logj_path, sep="\t")
    keep_lj = [
        "subject",
        "logJac_mean_MNI", "logJac_median_MNI", "logJac_p5_MNI", "logJac_p95_MNI",
        "logJac_mean_EBT", "logJac_median_EBT", "logJac_p5_EBT", "logJac_p95_EBT",
        "LesionVolume", "ScannerBrand", "ImageResolution"
    ]
    missing = [c for c in keep_lj if c not in lj.columns]
    if missing:
        raise RuntimeError(f"Missing columns in logJ TSV: {missing}\nHave: {list(lj.columns)}")
    lj = lj[keep_lj].copy()

    # lesion volume TSV (native): case column is the subject id
    les = pd.read_csv(lesion_path, sep="\t")
    if "case" not in les.columns or "Vlesion_native_mm3" not in les.columns:
        raise RuntimeError(f"lesion_fraction_qc_v2.tsv must contain case and Vlesion_native_mm3. Have: {list(les.columns)}")
    les = les.rename(columns={"case": "subject"})
    les = les[["subject", "Vlesion_native_mm3"]].copy()
    les["log1p_lesion_mm3"] = (les["Vlesion_native_mm3"].astype(float) + 1.0).apply(lambda x: __import__("math").log(x))

    # ventricle fraction TSV
    vf = pd.read_csv(vf_path, sep="\t")
    need_vf = ["subject", "ventricle_fraction", "brain_fraction"]
    missing = [c for c in need_vf if c not in vf.columns]
    if missing:
        raise RuntimeError(f"Missing columns in ventricle_brain_fraction.tsv: {missing}\nHave: {list(vf.columns)}")
    vf = vf[need_vf].copy()

    # merge
    df = cc.merge(lj, on="subject", how="inner") \
           .merge(les, on="subject", how="inner") \
           .merge(vf, on="subject", how="inner")

    # derived deltas (for convenience)
    df["dCC_harmonized"] = df["CC_EBT_bg"] - df["CC_MNI_bg"]
    if "gCC_MNI_bg" in df.columns and "gCC_EBT_bg" in df.columns:
        df["dgCC_harmonized"] = df["gCC_EBT_bg"] - df["gCC_MNI_bg"]
    df["dlogJac_mean"] = df["logJac_mean_EBT"] - df["logJac_mean_MNI"]

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print("[INFO] wrote:", out_path)
    print("[INFO] n:", df.shape[0], "cols:", df.shape[1])

if __name__ == "__main__":
    main()