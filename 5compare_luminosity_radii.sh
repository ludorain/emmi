#!/usr/bin/env bash

# ============================================================
# COMPARE LUMINOSITY RESULTS FOR R = 16, 20, 24 px
#
# Usage:
#   ./5compare_luminosity_radii.sh SENSOR CONSTANT PHASE
#
# Examples:
#   ./5compare_luminosity_radii.sh A1 T before_annealing
#   ./5compare_luminosity_radii.sh A1 T annealing_T=75_h=5
#   ./5compare_luminosity_radii.sh B1 v before_annealing
#
# CONSTANT:
#   T -> T is fixed, therefore luminosity is plotted vs overvoltage
#   v -> v is fixed, therefore luminosity is plotted vs temperature
#
# Nominal dataset:
#   R = 20 px
#   - statistical uncertainties from the R=20 master CSV
#   - luminosity systematic uncertainty deltaL already written in-place
#     by the systematic-uncertainty pipeline
#   - B/lambda systematic uncertainty read from spot_luminosity_final/anal
#
# Comparison datasets:
#   R = 16 px and R = 24 px
#   - experimental points with statistical errors
#   - independent fit, drawn as dashed line
#
# Rbest, R=15 and R=25 are intentionally NOT used.
# ============================================================

set -Eeuo pipefail
shopt -s nullglob

# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

R16_DATA_DIR="${R16_DATA_DIR:-$BASE_DIR/DATA_irradiated_isolated_R=16}"
R20_DATA_DIR="${R20_DATA_DIR:-$BASE_DIR/DATA_irradiated_isolated_R=20}"
R24_DATA_DIR="${R24_DATA_DIR:-$BASE_DIR/DATA_irradiated_isolated_R=24}"

SPOT_LUM_DIR="${SPOT_LUM_DIR:-$BASE_DIR/spot_luminosity_final}"
COMPARISON_DIR="${COMPARISON_DIR:-$BASE_DIR/comparison_results_radii}"

MACRO_V="$SPOT_LUM_DIR/lum_vs_v_fit.C"
MACRO_T="$SPOT_LUM_DIR/lum_vs_T_fit.C"

# ------------------------------------------------------------
# LOGGING / ERRORS
# ------------------------------------------------------------

log() {
    printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

on_error() {
    local exit_code=$?
    local line_no=$1
    printf '\nERROR: radius-comparison pipeline stopped at line %s (exit code %s).\n' \
        "$line_no" "$exit_code" >&2
    exit "$exit_code"
}

trap 'on_error $LINENO' ERR

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------

if [[ $# -ne 3 ]]; then
    cat >&2 <<EOF
Usage:
  $0 SENSOR CONSTANT PHASE

Examples:
  $0 A1 T before_annealing
  $0 A1 T annealing_T=75_h=5
  $0 B1 v before_annealing
EOF
    exit 1
fi

SENSOR="$1"
CONSTANT="$2"
PHASE="$3"

case "$SENSOR" in
    A1|A2|B1|B2) ;;
    *) die "Invalid sensor '$SENSOR'. Allowed values: A1 A2 B1 B2." ;;
esac

case "$CONSTANT" in
    T|v) ;;
    *) die "Invalid constant '$CONSTANT'. Allowed values: T or v." ;;
esac

# ------------------------------------------------------------
# REQUIREMENTS
# ------------------------------------------------------------

command -v root >/dev/null 2>&1 || die "ROOT executable 'root' not found in PATH."

for d in "$R16_DATA_DIR" "$R20_DATA_DIR" "$R24_DATA_DIR"; do
    [[ -d "$d" ]] || die "Required data directory not found: $d"
done

[[ -d "$SPOT_LUM_DIR" ]] || die "Directory not found: $SPOT_LUM_DIR"
[[ -f "$MACRO_V" ]] || die "ROOT macro not found: $MACRO_V"
[[ -f "$MACRO_T" ]] || die "ROOT macro not found: $MACRO_T"

mkdir -p "$COMPARISON_DIR"

# ------------------------------------------------------------
# DISCOVER THE FIXED VALUE FROM THE R=20 MERGED DIRECTORY
#
# Example:
#   SENSOR=A1, CONSTANT=T  -> A1_T=20
#   SENSOR=A1, CONSTANT=v  -> A1_v=5
# ------------------------------------------------------------

MERGED20_ROOT="$R20_DATA_DIR/merged_files"
[[ -d "$MERGED20_ROOT" ]] || die "Merged-files directory not found: $MERGED20_ROOT"

RUN_MATCHES=()
while IFS= read -r -d '' d; do
    RUN_MATCHES+=("$d")
done < <(
    find "$MERGED20_ROOT" \
        -mindepth 1 -maxdepth 1 \
        -type d \
        -name "${SENSOR}_${CONSTANT}=*" \
        -print0
)

if (( ${#RUN_MATCHES[@]} == 0 )); then
    die "No merged directory found for ${SENSOR}_${CONSTANT}=* in $MERGED20_ROOT"
fi

if (( ${#RUN_MATCHES[@]} > 1 )); then
    printf 'Matching merged directories:\n' >&2
    printf '  %s\n' "${RUN_MATCHES[@]}" >&2
    die "Expected exactly one fixed-value directory for ${SENSOR}_${CONSTANT}."
fi

RUN_PREFIX="$(basename "${RUN_MATCHES[0]}")"
FIXED_VALUE="${RUN_PREFIX#${SENSOR}_${CONSTANT}=}"

# ------------------------------------------------------------
# MASTER CSV FILES
# ------------------------------------------------------------

MASTER_NAME="${SENSOR}_${CONSTANT}_all_phases_global_ID.csv"

R16_MERGED_DIR="$R16_DATA_DIR/merged_files/$RUN_PREFIX"
R20_MERGED_DIR="$R20_DATA_DIR/merged_files/$RUN_PREFIX"
R24_MERGED_DIR="$R24_DATA_DIR/merged_files/$RUN_PREFIX"

R16_CSV="$R16_MERGED_DIR/$MASTER_NAME"
R20_CSV="$R20_MERGED_DIR/$MASTER_NAME"
R24_CSV="$R24_MERGED_DIR/$MASTER_NAME"

for f in "$R16_CSV" "$R20_CSV" "$R24_CSV"; do
    [[ -f "$f" ]] || die "Required all-phases CSV not found: $f"
done

# The R=20 file must already contain the point-by-point systematic uncertainty
# produced by 4luminosity_sistematical_uncertainty.sh.
R20_HEADER="$(head -n 1 "$R20_CSV" | tr -d '\r')"

if [[ ",${R20_HEADER}," != *",deltaL,"* ]]; then
    die "The nominal R=20 master CSV does not contain the 'deltaL' column: $R20_CSV"
fi

# ------------------------------------------------------------
# SYSTEMATIC FIT-PARAMETER FILE FROM 5analysis
#
# The exact prefix can evolve, therefore the pipeline looks inside the unique
# phase output directory and requires exactly one *_B_values.csv or
# *_lambda_values.csv file.
# ------------------------------------------------------------

ANALYSIS_PHASE_DIR="$SPOT_LUM_DIR/analysis/$RUN_PREFIX/$PHASE"
[[ -d "$ANALYSIS_PHASE_DIR" ]] || die \
    "5analysis output directory not found: $ANALYSIS_PHASE_DIR"

if [[ "$CONSTANT" == "T" ]]; then
    SYS_PATTERN="*_B_values.csv"
    PARAMETER_NAME="B"
else
    SYS_PATTERN="*_lambda_values.csv"
    PARAMETER_NAME="lambda"
fi

SYS_MATCHES=( "$ANALYSIS_PHASE_DIR"/$SYS_PATTERN )

if (( ${#SYS_MATCHES[@]} == 0 )); then
    die "No $PARAMETER_NAME systematic-values CSV matching '$SYS_PATTERN' found in $ANALYSIS_PHASE_DIR"
fi

if (( ${#SYS_MATCHES[@]} > 1 )); then
    printf 'Matching systematic-value files:\n' >&2
    printf '  %s\n' "${SYS_MATCHES[@]}" >&2
    die "Expected exactly one $PARAMETER_NAME systematic-values CSV in $ANALYSIS_PHASE_DIR"
fi

SYSTEMATIC_VALUES_CSV="${SYS_MATCHES[0]}"

# ------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------

OUTPUT_DIR="$COMPARISON_DIR/$RUN_PREFIX/$PHASE"
mkdir -p "$OUTPUT_DIR"

# Remove only the single-hotspot comparison products, so stale plots from an
# older radius comparison cannot survive a rerun.
rm -rf \
    "$OUTPUT_DIR/lum_vs_v_all_spots" \
    "$OUTPUT_DIR/lum_vs_T_all_spots_expfit"

PREFIX="${RUN_PREFIX}_${PHASE}_R16_R20_R24"

# ------------------------------------------------------------
# REPORT
# ------------------------------------------------------------

printf '\n============================================================\n'
printf 'LUMINOSITY COMPARISON: R = 16, 20, 24 px\n'
printf 'Sensor:          %s\n' "$SENSOR"
printf 'Fixed condition: %s=%s\n' "$CONSTANT" "$FIXED_VALUE"
printf 'Phase:           %s\n' "$PHASE"
printf 'Nominal radius:  R = 20 px\n'
printf 'Output:          %s\n' "$OUTPUT_DIR"
printf '============================================================\n'

printf '\nInput files:\n'
printf '  R = 16: %s\n' "$R16_CSV"
printf '  R = 20: %s\n' "$R20_CSV"
printf '  R = 24: %s\n' "$R24_CSV"
printf '  %s systematics: %s\n' "$PARAMETER_NAME" "$SYSTEMATIC_VALUES_CSV"

# ------------------------------------------------------------
# ROOT MACRO CALL
# ------------------------------------------------------------

if [[ "$CONSTANT" == "T" ]]; then
    log "Drawing luminosity vs overvoltage with R=16/20/24"

    root -l -b -q \
        "${MACRO_V}(\"${R20_CSV}\",\"${PHASE}\",\"${SYSTEMATIC_VALUES_CSV}\",\"${OUTPUT_DIR}\",\"${PREFIX}\",\"${R16_CSV}\",\"${R24_CSV}\")"
else
    log "Drawing luminosity vs temperature with R=16/20/24"

    root -l -b -q \
        "${MACRO_T}(\"${R20_CSV}\",\"${PHASE}\",\"${SYSTEMATIC_VALUES_CSV}\",\"${OUTPUT_DIR}\",\"${PREFIX}\",\"${R16_CSV}\",\"${R24_CSV}\")"
fi

printf '\n============================================================\n'
printf 'RADIUS COMPARISON COMPLETED SUCCESSFULLY\n'
printf 'Plots saved in:\n  %s\n' "$OUTPUT_DIR"
printf '============================================================\n'
