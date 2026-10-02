#!/usr/bin/env bash

set -u
set -o pipefail

# ============================================================
# Run lum_vs_T_fit.C for ONE sensor and ONE annealing phase.
#
# Usage:
#   ./run_lum_vs_T_single_phase.sh SENSOR PHASE
#
# Examples:
#   ./run_lum_vs_T_single_phase.sh A1 before_annealing
#   ./run_lum_vs_T_single_phase.sh A1 annealing_T=75_h=5
#   ./run_lum_vs_T_single_phase.sh B1 annealing_T=150_h=25
#
# ROOT is intentionally left open after the macro execution.
# ============================================================

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 SENSOR PHASE"
    exit 1
fi

SENSOR="$1"
PHASE="$2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MACRO_DIR="$SCRIPT_DIR/spot_luminosity_final"
R20_DIR="$SCRIPT_DIR/DATA_irradiated_isolated_R=20"


# ============================================================
# Fixed-voltage condition
# ============================================================

case "$SENSOR" in
    B1)
        COND="${SENSOR}_v=3"
        ;;
    A1|A2|B2)
        COND="${SENSOR}_v=5"
        ;;
    *)
        echo "ERROR: invalid sensor '$SENSOR'."
        exit 2
        ;;
esac


# ============================================================
# Input / output paths
# ============================================================

ALL20="$R20_DIR/merged_files/$COND/${SENSOR}_v_all_phases_global_ID.csv"

OUTDIR="$MACRO_DIR/analysis/$COND/$PHASE"
PREFIX="${COND}_${PHASE}"

SYSCsv="$OUTDIR/${PREFIX}_lambda_values.csv"

SUMMARY="$MACRO_DIR/analysis/$COND/${COND}_lambda_vs_phase.csv"


# ============================================================
# Checks
# ============================================================

if [[ ! -f "$ALL20" ]]; then
    echo "ERROR: input file not found:"
    echo "$ALL20"
    exit 3
fi

mkdir -p "$OUTDIR"


# If the systematic-lambda file does not already exist,
# run the nominal fit without it.
if [[ ! -f "$SYSCsv" ]]; then
    echo "WARNING: systematic file not found:"
    echo "  $SYSCsv"
    echo "Running lum_vs_T_fit.C without lambda systematic uncertainties."
    SYSCsv=""
fi


# ============================================================
# Run ROOT
#
# No -q  -> ROOT remains open
# No -b  -> graphical canvases remain available
# ============================================================

cd "$MACRO_DIR" || exit 1

echo
echo "Sensor : $SENSOR"
echo "Phase  : $PHASE"
echo "Input  : $ALL20"
echo

root -l \
"lum_vs_T_fit.C(\"$ALL20\",\"$PHASE\",\"$SYSCsv\",\"$OUTDIR\",\"$PREFIX\",\"\",\"\",\"$SUMMARY\")"