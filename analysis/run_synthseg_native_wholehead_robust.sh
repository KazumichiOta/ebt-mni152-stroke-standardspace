#!/usr/bin/env bash
#
# run_synthseg_native_wholehead_robust.sh
#
# Purpose:
#   Generate native-space SynthSeg labels for all ATLAS subjects using
#   the N4-corrected, non-skull-stripped whole-head T1 image.
#
#       anat/T1w_n4.nii.gz
#
#   Outputs are written to a NEW directory:
#
#       anat/synthseg_native_wholehead_robust/labels.nii.gz
#
#   SynthSeg is run with:
#
#       --parc --robust
#
# Environment variables:
#
#   ROOT=/path/to/ATLAS_workdir
#   SYNTHSEG_CMD=mri_synthseg
#   JOBS=2
#   THREADS=<automatically calculated>
#   OVERWRITE=0
#   LIMIT=0
#
# LIMIT=0 means all subjects.
# Example test:
#
#   ROOT=/path/to/ATLAS_workdir LIMIT=3 JOBS=1 \
#       bash analysis/run_synthseg_native_wholehead_robust.sh
#

set -euo pipefail

export LC_ALL=C
export LANG=C


# =========================================================================
# Configuration
# =========================================================================

ROOT="${ROOT:-}"

SYNTHSEG_CMD="${SYNTHSEG_CMD:-mri_synthseg}"

INPUT_NAME="T1w_n4.nii.gz"

OUTPUT_DIR_NAME="synthseg_native_wholehead_robust"
OUTPUT_NAME="labels.nii.gz"

# Number of subjects processed simultaneously.
JOBS="${JOBS:-2}"

# 0 = skip existing new outputs
# 1 = overwrite existing new outputs
OVERWRITE="${OVERWRITE:-0}"

# 0 = all subjects
# positive integer = only first N subjects, useful for testing
LIMIT="${LIMIT:-0}"


# =========================================================================
# CPU configuration
# =========================================================================

HW_CPUS="$(/usr/sbin/sysctl -n hw.ncpu 2>/dev/null || echo 4)"

if [[ -z "${THREADS:-}" ]]; then
    THREADS=$(( HW_CPUS / JOBS ))

    if [[ "$THREADS" -lt 1 ]]; then
        THREADS=1
    fi
fi


# =========================================================================
# Checks
# =========================================================================

if [[ -z "$ROOT" ]]; then
    echo "[FATAL] ROOT is not set."
    echo "        Example:"
    echo "        ROOT=/path/to/ATLAS_workdir bash $0"
    exit 1
fi

if [[ ! -d "$ROOT" ]]; then
    echo "[FATAL] ROOT does not exist:"
    echo "        $ROOT"
    exit 1
fi

if [[ "$SYNTHSEG_CMD" == */* ]]; then

    if [[ ! -x "$SYNTHSEG_CMD" ]]; then
        echo "[FATAL] mri_synthseg is not executable:"
        echo "        $SYNTHSEG_CMD"
        exit 1
    fi

else

    SYNTHSEG_RESOLVED="$(command -v "$SYNTHSEG_CMD" 2>/dev/null || true)"

    if [[ -z "$SYNTHSEG_RESOLVED" ]]; then
        echo "[FATAL] mri_synthseg was not found in PATH."
        echo "        Set SYNTHSEG_CMD explicitly if needed."
        exit 1
    fi

    SYNTHSEG_CMD="$SYNTHSEG_RESOLVED"

fi

if [[ "$JOBS" -lt 1 ]]; then
    echo "[FATAL] JOBS must be >= 1"
    exit 1
fi

if [[ "$THREADS" -lt 1 ]]; then
    echo "[FATAL] THREADS must be >= 1"
    exit 1
fi


# =========================================================================
# Logging
# =========================================================================

QC_DIR="${ROOT}/_qc/synthseg_native_wholehead_robust"

mkdir -p "$QC_DIR"

RUN_ID="$(date '+%Y%m%d_%H%M%S')"

STATUS_FILE="${QC_DIR}/status_${RUN_ID}.tsv"
SUBJECT_LIST="${QC_DIR}/subjects_${RUN_ID}.txt"

echo -e "status\tsubject\tinput\toutput\tnote" > "$STATUS_FILE"


# =========================================================================
# Build subject list
# =========================================================================

: > "$SUBJECT_LIST"

for SUB_DIR in "$ROOT"/sub-*; do

    [[ -d "$SUB_DIR" ]] || continue

    SID="$(basename "$SUB_DIR")"

    echo "$SID" >> "$SUBJECT_LIST"

done

sort -o "$SUBJECT_LIST" "$SUBJECT_LIST"


N_FOUND="$(wc -l < "$SUBJECT_LIST" | tr -d ' ')"

if [[ "$N_FOUND" -eq 0 ]]; then
    echo "[FATAL] No sub-* directories found."
    exit 1
fi


if [[ "$LIMIT" -gt 0 ]]; then

    SUBJECT_LIST_LIMITED="${QC_DIR}/subjects_${RUN_ID}_limited.txt"

    head -n "$LIMIT" "$SUBJECT_LIST" > "$SUBJECT_LIST_LIMITED"

    SUBJECT_LIST="$SUBJECT_LIST_LIMITED"

fi


N_RUN="$(wc -l < "$SUBJECT_LIST" | tr -d ' ')"


# =========================================================================
# Run information
# =========================================================================

echo
echo "============================================================"
echo "SynthSeg robust whole-head native-space segmentation"
echo "============================================================"
echo
echo "[INFO] ROOT:"
echo "       $ROOT"
echo
echo "[INFO] INPUT:"
echo "       anat/$INPUT_NAME"
echo
echo "[INFO] OUTPUT:"
echo "       anat/$OUTPUT_DIR_NAME/$OUTPUT_NAME"
echo
echo "[INFO] SynthSeg:"
echo "       $SYNTHSEG_CMD"
echo
echo "[INFO] Subjects found : $N_FOUND"
echo "[INFO] Subjects queued: $N_RUN"
echo "[INFO] JOBS           : $JOBS"
echo "[INFO] THREADS/job    : $THREADS"
echo "[INFO] Hardware CPUs  : $HW_CPUS"
echo "[INFO] OVERWRITE      : $OVERWRITE"
echo "[INFO] LIMIT          : $LIMIT"
echo
echo "[INFO] Status file:"
echo "       $STATUS_FILE"
echo


# =========================================================================
# Process one subject
# =========================================================================

run_one() {

    SID="$1"

    ANAT_DIR="${ROOT}/${SID}/anat"

    IN_IMG="${ANAT_DIR}/${INPUT_NAME}"

    OUT_DIR="${ANAT_DIR}/${OUTPUT_DIR_NAME}"

    OUT_LBL="${OUT_DIR}/${OUTPUT_NAME}"

    LOG_FILE="${QC_DIR}/${RUN_ID}_${SID}.log"


    # ---------------------------------------------------------------------
    # Missing input
    # ---------------------------------------------------------------------

    if [[ ! -f "$IN_IMG" ]]; then

        echo "[SKIP] $SID : missing $INPUT_NAME"

        printf "MISSING_INPUT\t%s\t%s\t%s\t%s\n" \
            "$SID" \
            "$IN_IMG" \
            "$OUT_LBL" \
            "input_not_found" \
            >> "$STATUS_FILE"

        return 0
    fi


    # ---------------------------------------------------------------------
    # Existing NEW output
    # ---------------------------------------------------------------------

    if [[ -f "$OUT_LBL" && "$OVERWRITE" != "1" ]]; then

        echo "[SKIP] $SID : whole-head SynthSeg output already exists"

        printf "EXISTS\t%s\t%s\t%s\t%s\n" \
            "$SID" \
            "$IN_IMG" \
            "$OUT_LBL" \
            "not_overwritten" \
            >> "$STATUS_FILE"

        return 0
    fi


    mkdir -p "$OUT_DIR"


    if [[ "$OVERWRITE" == "1" && -f "$OUT_LBL" ]]; then
        rm -f "$OUT_LBL"
    fi


    # ---------------------------------------------------------------------
    # SynthSeg
    #
    # Important:
    #   No HD-BET
    #   No 0-1 normalization
    #   No BM4D
    #
    # Input is the N4-corrected whole-head native T1.
    #
    # Final analysis pipeline:
    #   N4-corrected whole-head native T1
    #   -> SynthSeg --parc --robust.
    # ---------------------------------------------------------------------

    echo "[RUN ] $SID"

    {
        echo "============================================================"
        echo "Subject: $SID"
        echo "Input  : $IN_IMG"
        echo "Output : $OUT_LBL"
        echo "Command:"
        echo "\"$SYNTHSEG_CMD\" --i \"$IN_IMG\" --o \"$OUT_LBL\" --parc --robust --threads \"$THREADS\""
        echo "============================================================"
        echo
    } > "$LOG_FILE"


  if "$SYNTHSEG_CMD" \
    --i "$IN_IMG" \
    --o "$OUT_LBL" \
    --parc \
    --robust \
    --threads "$THREADS" \
    >> "$LOG_FILE" 2>&1
then
        if [[ -s "$OUT_LBL" ]]; then

            echo "[ OK ] $SID"

            printf "OK\t%s\t%s\t%s\t%s\n" \
                "$SID" \
                "$IN_IMG" \
                "$OUT_LBL" \
                "completed" \
                >> "$STATUS_FILE"

        else

            echo "[ERROR] $SID : command completed but output is missing/empty"

            printf "ERROR\t%s\t%s\t%s\t%s\n" \
                "$SID" \
                "$IN_IMG" \
                "$OUT_LBL" \
                "output_missing_or_empty" \
                >> "$STATUS_FILE"
        fi

    else

        echo "[ERROR] $SID : mri_synthseg failed"

        printf "ERROR\t%s\t%s\t%s\t%s\n" \
            "$SID" \
            "$IN_IMG" \
            "$OUT_LBL" \
            "mri_synthseg_failed" \
            >> "$STATUS_FILE"
    fi
}


export -f run_one

export ROOT
export SYNTHSEG_CMD
export INPUT_NAME
export OUTPUT_DIR_NAME
export OUTPUT_NAME
export OVERWRITE
export THREADS
export QC_DIR
export RUN_ID
export STATUS_FILE


# =========================================================================
# Parallel execution
# =========================================================================

cat "$SUBJECT_LIST" \
    | xargs -P "$JOBS" -I{} bash -c 'run_one "$1"' _ {}


# =========================================================================
# Summary
# =========================================================================

N_OK="$(awk -F'\t' '$1=="OK"{n++} END{print n+0}' "$STATUS_FILE")"

N_EXISTS="$(awk -F'\t' '$1=="EXISTS"{n++} END{print n+0}' "$STATUS_FILE")"

N_MISSING="$(awk -F'\t' '$1=="MISSING_INPUT"{n++} END{print n+0}' "$STATUS_FILE")"

N_ERROR="$(awk -F'\t' '$1=="ERROR"{n++} END{print n+0}' "$STATUS_FILE")"


echo
echo "============================================================"
echo "Finished"
echo "============================================================"
echo
echo "OK             : $N_OK"
echo "Already exists : $N_EXISTS"
echo "Missing input  : $N_MISSING"
echo "Errors         : $N_ERROR"
echo "Queued         : $N_RUN"
echo
echo "Status:"
echo "$STATUS_FILE"
echo
echo "Outputs:"
echo "$ROOT/sub-*/anat/$OUTPUT_DIR_NAME/$OUTPUT_NAME"
echo


if [[ "$N_ERROR" -gt 0 ]]; then
    echo "[WARN] One or more SynthSeg jobs failed."
    echo "       Check individual logs under:"
    echo "       $QC_DIR"
    exit 2
fi

echo "[DONE] SynthSeg robust whole-head segmentation completed."

