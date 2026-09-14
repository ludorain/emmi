#!/usr/bin/env bash

# ============================================================
# COMPARE LUMINOSITY RESULTS FOR DIFFERENT INTEGRATION RADII
#
# Usage:
#   ./compare_luminosity_radii.sh SENSOR CONSTANT PHASE [GLOBAL_SPOT_ID]
#
# Examples:
#   ./compare_luminosity_radii.sh A1 T before_annealing 31
#   ./compare_luminosity_radii.sh A1 T annealing_T=75_h=5
#   ./compare_luminosity_radii.sh B1 v before_annealing 12
#
# Arguments:
#   SENSOR          A1 | A2 | B1 | B2
#   CONSTANT        T | v
#   PHASE           before_annealing or any annealing_* directory
#   GLOBAL_SPOT_ID  optional; if omitted, all global hotspot IDs are plotted
#
# The script must be located inside the main emmi directory.
# ROOT macros must be in:
#   emmi/comparison_results_radii/
#       lum_vs_V_comparison.C
#       lum_vs_T_comparison.C
#
# Compared datasets:
#   R best  -> DATA_irradiated_isolated_changeR
#   R = 15  -> DATA_irradiated_isolated_R=15
#   R = 20  -> DATA_irradiated_isolated_R=20
#   R = 25  -> DATA_irradiated_isolated_R=25
# ============================================================

set -Eeuo pipefail
shopt -s nullglob

# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BEST_DATA_DIR="${BEST_DATA_DIR:-$BASE_DIR/DATA_irradiated_isolated_changeR}"
R15_DATA_DIR="${R15_DATA_DIR:-$BASE_DIR/DATA_irradiated_isolated_R=15}"
R20_DATA_DIR="${R20_DATA_DIR:-$BASE_DIR/DATA_irradiated_isolated_R=20}"
R25_DATA_DIR="${R25_DATA_DIR:-$BASE_DIR/DATA_irradiated_isolated_R=25}"

COMPARISON_DIR="${COMPARISON_DIR:-$BASE_DIR/comparison_results_radii}"
MACRO_V="$COMPARISON_DIR/lum_vs_V_comparison.C"
MACRO_T="$COMPARISON_DIR/lum_vs_T_comparison.C"

# ------------------------------------------------------------
# LOGGING / ERRORS
# ------------------------------------------------------------

log() {
    printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"
}

warn() {
    printf 'WARNING: %s\n' "$*" >&2
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

on_error() {
    local exit_code=$?
    local line_no=$1
    printf '\nERROR: comparison pipeline stopped at line %s (exit code %s).\n' \
        "$line_no" "$exit_code" >&2
    exit "$exit_code"
}

trap 'on_error $LINENO' ERR

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------

if [[ $# -ne 3 && $# -ne 4 ]]; then
    cat >&2 <<EOF
Usage:
  $0 SENSOR CONSTANT PHASE [GLOBAL_SPOT_ID]

Examples:
  $0 A1 T before_annealing 31
  $0 A1 T annealing_T=75_h=5
  $0 B1 v before_annealing 12
EOF
    exit 1
fi

SENSOR="$1"
CONSTANT="$2"
PHASE="$3"
SELECTED_SPOT="${4:--1}"

case "$SENSOR" in
    A1|A2|B1|B2) ;;
    *) die "Invalid sensor '$SENSOR'. Allowed values: A1 A2 B1 B2." ;;
esac

case "$CONSTANT" in
    T|v) ;;
    *) die "Invalid constant '$CONSTANT'. Allowed values: T or v." ;;
esac

if [[ "$SELECTED_SPOT" != "-1" && ! "$SELECTED_SPOT" =~ ^[0-9]+$ ]]; then
    die "GLOBAL_SPOT_ID must be a non-negative integer. Received: $SELECTED_SPOT"
fi

# ------------------------------------------------------------
# REQUIREMENTS
# ------------------------------------------------------------

command -v root >/dev/null 2>&1 || die "ROOT executable 'root' not found in PATH."

for d in "$BEST_DATA_DIR" "$R15_DATA_DIR" "$R20_DATA_DIR" "$R25_DATA_DIR"; do
    [[ -d "$d" ]] || die "Required data directory not found: $d"
done

mkdir -p "$COMPARISON_DIR"

if [[ "$CONSTANT" == "T" ]]; then
    [[ -f "$MACRO_V" ]] || die "ROOT macro not found: $MACRO_V"
else
    [[ -f "$MACRO_T" ]] || die "ROOT macro not found: $MACRO_T"
fi

# ------------------------------------------------------------
# DISCOVER THE RUN FROM THE BEST-RADIUS REFERENCE DATASET
# ------------------------------------------------------------

SOURCE_PHASE_DIR="$BEST_DATA_DIR/$PHASE"
[[ -d "$SOURCE_PHASE_DIR" ]] || die "Phase directory not found: $SOURCE_PHASE_DIR"

RUN_MATCHES=()
while IFS= read -r -d '' d; do
    RUN_MATCHES+=("$d")
done < <(
    find "$SOURCE_PHASE_DIR" \
        -mindepth 1 -maxdepth 1 \
        -type d \
        -name "${SENSOR}_${CONSTANT}=*_run=*" \
        -print0
)

if (( ${#RUN_MATCHES[@]} == 0 )); then
    die "No run found for ${SENSOR}_${CONSTANT}=* in phase '$PHASE'."
fi

if (( ${#RUN_MATCHES[@]} > 1 )); then
    printf 'Matching runs found:\n' >&2
    printf '  %s\n' "${RUN_MATCHES[@]}" >&2
    die "Expected exactly one run for the selected sensor/constant/phase."
fi

SOURCE_RUN_DIR="${RUN_MATCHES[0]}"
RUN_NAME="$(basename "$SOURCE_RUN_DIR")"

if [[ ! "$RUN_NAME" =~ ^(${SENSOR}_${CONSTANT}=([^_]+))_run=(.+)$ ]]; then
    die "Cannot parse run directory name: $RUN_NAME"
fi

RUN_PREFIX="${BASH_REMATCH[1]}"
FIXED_VALUE="${BASH_REMATCH[2]}"
RUN_NUMBER="${BASH_REMATCH[3]}"
FIXED_LABEL="${CONSTANT}=${FIXED_VALUE}"

# ------------------------------------------------------------
# INPUT CSV PATHS
# ------------------------------------------------------------

BEST_MERGED_DIR="$BEST_DATA_DIR/merged_files/$RUN_PREFIX"
R15_MERGED_DIR="$R15_DATA_DIR/merged_files/$RUN_PREFIX"
R20_MERGED_DIR="$R20_DATA_DIR/merged_files/$RUN_PREFIX"
R25_MERGED_DIR="$R25_DATA_DIR/merged_files/$RUN_PREFIX"

FILE_STEM="${RUN_PREFIX}_${PHASE}_run=${RUN_NUMBER}"

# Prefer the global_complete CSV for the best-radius dataset as well,
# because it contains the integration_radius information needed to plot
# the fit parameter versus radius. Fall back to the historical global_ID
# file only if global_complete is not available.
BEST_CSV="$BEST_MERGED_DIR/${FILE_STEM}_global_complete.csv"
if [[ ! -f "$BEST_CSV" ]]; then
    BEST_CSV="$BEST_MERGED_DIR/${FILE_STEM}_global_ID.csv"
fi

R15_CSV="$R15_MERGED_DIR/${FILE_STEM}_global_complete.csv"
R20_CSV="$R20_MERGED_DIR/${FILE_STEM}_global_complete.csv"
R25_CSV="$R25_MERGED_DIR/${FILE_STEM}_global_complete.csv"

for f in "$BEST_CSV" "$R15_CSV" "$R20_CSV" "$R25_CSV"; do
    [[ -f "$f" ]] || die "Required comparison CSV not found: $f"
done

# ------------------------------------------------------------
# OUTPUT DIRECTORY
# ------------------------------------------------------------

OUTPUT_DIR="$COMPARISON_DIR/$PHASE/$RUN_NAME"
mkdir -p "$OUTPUT_DIR"

# Cleaning policy:
# - if no specific spot is requested, remove all previously generated plots
#   in the output directory and regenerate the full set;
# - if a specific spot is requested, remove only files associated with that
#   spot ID, leaving the plots of other spots untouched.
if [[ "$SELECTED_SPOT" == "-1" ]]; then
    find "$OUTPUT_DIR" -maxdepth 1 -type f         \( -name '*.png' -o -name '*.pdf' \)         -delete
else
    find "$OUTPUT_DIR" -maxdepth 1 -type f         \( -name "*spot${SELECTED_SPOT}*.png" -o -name "*spot${SELECTED_SPOT}*.pdf" \)         -delete
fi

# ------------------------------------------------------------
# REPORT
# ------------------------------------------------------------

printf '\n============================================================\n'
printf 'LUMINOSITY COMPARISON BETWEEN INTEGRATION RADII\n'
printf 'Sensor:          %s\n' "$SENSOR"
printf 'Constant:        %s=%s\n' "$CONSTANT" "$FIXED_VALUE"
printf 'Phase:           %s\n' "$PHASE"
printf 'Run:             %s\n' "$RUN_NUMBER"
if [[ "$SELECTED_SPOT" == "-1" ]]; then
    printf 'Global spot ID:  ALL\n'
else
    printf 'Global spot ID:  %s\n' "$SELECTED_SPOT"
fi
printf 'Output:          %s\n' "$OUTPUT_DIR"
printf '============================================================\n'

printf '\nInput files:\n'
printf '  R best: %s\n' "$BEST_CSV"
printf '  R = 15: %s\n' "$R15_CSV"
printf '  R = 20: %s\n' "$R20_CSV"
printf '  R = 25: %s\n' "$R25_CSV"

# ------------------------------------------------------------
# ROOT MACRO CALL
# ------------------------------------------------------------

if [[ "$CONSTANT" == "T" ]]; then
    log "Running luminosity vs overvoltage comparison"

    root -l -b -q \
        "${MACRO_V}(\"${BEST_CSV}\",\"${R15_CSV}\",\"${R20_CSV}\",\"${R25_CSV}\",\"${OUTPUT_DIR}\",\"${SENSOR}\",\"${FIXED_LABEL}\",\"${PHASE}\",${SELECTED_SPOT})"
else
    log "Running luminosity vs temperature comparison"

    root -l -b -q \
        "${MACRO_T}(\"${BEST_CSV}\",\"${R15_CSV}\",\"${R20_CSV}\",\"${R25_CSV}\",\"${OUTPUT_DIR}\",\"${SENSOR}\",\"${FIXED_LABEL}\",\"${PHASE}\",${SELECTED_SPOT})"
fi

printf '\n============================================================\n'
printf 'COMPARISON COMPLETED SUCCESSFULLY\n'
printf 'Plots saved in:\n  %s\n' "$OUTPUT_DIR"
printf '============================================================\n'
