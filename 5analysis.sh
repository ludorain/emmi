#!/usr/bin/env bash
set -u
set -o pipefail

# ============================================================
# 5analysis.sh
# Final EMMI analysis driver.
#
# Usage:
#   ./5analysis.sh SENSOR [T|v] [analysis|phases]
#   ./5analysis.sh SENSOR [analysis|phases]   # both constants
#   ./5analysis.sh SENSOR [T|v]              # both commands
#   ./5analysis.sh SENSOR                    # both constants and commands
#
# SENSOR is mandatory. A call with no arguments is rejected intentionally.
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO_DIR="$SCRIPT_DIR/spot_luminosity_final"

R16_DIR="$SCRIPT_DIR/DATA_irradiated_isolated_R=16"
R20_DIR="$SCRIPT_DIR/DATA_irradiated_isolated_R=20"
R24_DIR="$SCRIPT_DIR/DATA_irradiated_isolated_R=24"

PHASES=(
  "before_annealing"
  "annealing_T=75_h=5"
  "annealing_T=75_h=25"
  "annealing_T=100_h=5"
  "annealing_T=100_h=25"
  "annealing_T=125_h=5"
  "annealing_T=125_h=25"
  "annealing_T=150_h=5"
  "annealing_T=150_h=25"
)

usage() {
  cat <<USAGE
Usage: $0 SENSOR [T|v] [analysis|phases]
       $0 SENSOR [analysis|phases]

Examples:
  $0 A1 T analysis
  $0 B2 phases
  $0 A2 v
  $0 A1
USAGE
}

if [[ $# -eq 0 ]]; then
  echo "ERROR: at least the sensor must be specified." >&2
  usage >&2
  exit 2
fi

SENSOR="$1"
shift
case "$SENSOR" in
  A1|A2|B1|B2) ;;
  *) echo "ERROR: invalid sensor '$SENSOR'. Expected A1, A2, B1 or B2." >&2; exit 2 ;;
esac

CONSTANT=""
COMMAND=""
for arg in "$@"; do
  case "$arg" in
    T|v)
      [[ -z "$CONSTANT" ]] || { echo "ERROR: constant specified more than once." >&2; exit 2; }
      CONSTANT="$arg"
      ;;
    analysis|phases)
      [[ -z "$COMMAND" ]] || { echo "ERROR: command specified more than once." >&2; exit 2; }
      COMMAND="$arg"
      ;;
    *)
      echo "ERROR: unknown argument '$arg'." >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$CONSTANT" ]]; then CONSTANTS=(T v); else CONSTANTS=("$CONSTANT"); fi
if [[ -z "$COMMAND" ]]; then COMMANDS=(analysis phases); else COMMANDS=("$COMMAND"); fi

for f in analysis_common.h phase_common.h lum_vs_v_systematics.C lum_vs_T_systematics.C lum_vs_v_fit.C lum_vs_T_fit.C lum_vs_phase_T_const.C lum_vs_phase_v_const.C; do
  if [[ ! -f "$MACRO_DIR/$f" ]]; then
    echo "ERROR: missing macro/helper $MACRO_DIR/$f" >&2
    exit 3
  fi
done

if ! command -v root >/dev/null 2>&1; then
  echo "ERROR: ROOT executable 'root' was not found in PATH." >&2
  exit 3
fi

condition_name() {
  local c="$1"
  if [[ "$c" == "T" ]]; then printf '%s_T=20' "$SENSOR"; else printf '%s_v=5' "$SENSOR"; fi
}

all_phases_file() {
  local c="$1" cond
  cond="$(condition_name "$c")"
  if [[ "$c" == "T" ]]; then
    printf '%s/merged_files/%s/%s_T_all_phases_global_ID.csv' "$R20_DIR" "$cond" "$SENSOR"
  else
    printf '%s/merged_files/%s/%s_v_all_phases_global_ID.csv' "$R20_DIR" "$cond" "$SENSOR"
  fi
}

# Return the newest run (lexicographically, matching YYYYMMDD-HHMMSS) for one phase.
# Missing phases are allowed and return an empty string.
find_phase_file() {
  local base="$1" c="$2" ph="$3" cond dir pattern
  cond="$(condition_name "$c")"
  dir="$base/merged_files/$cond"
  pattern="$dir/${cond}_${ph}_run=*_global_complete.csv"
  compgen -G "$pattern" | sort | tail -n 1 || true
}

run_root() {
  local expr="$1"
  ( cd "$MACRO_DIR" && root -l -b -q "$expr" )
}

run_analysis() {
  local c="$1" cond all20 analysis_base
  cond="$(condition_name "$c")"
  all20="$(all_phases_file "$c")"
  analysis_base="$MACRO_DIR/analysis/$cond"

  if [[ ! -f "$all20" ]]; then
    echo "WARNING: nominal R=20 all-phases file not found: $all20" >&2
    return 0
  fi

  for ph in "${PHASES[@]}"; do
    local f16 f24 outdir prefix syscsv
    f16="$(find_phase_file "$R16_DIR" "$c" "$ph")"
    f24="$(find_phase_file "$R24_DIR" "$c" "$ph")"

    # A phase absent from the nominal all-phases CSV will be skipped by the ROOT macro.
    outdir="$analysis_base/$ph"
    prefix="${cond}_${ph}"

    # Recreate the complete sensor/constant/phase output to prevent stale files.
    rm -rf "$outdir"
    mkdir -p "$outdir"

    if [[ "$c" == "T" ]]; then
      syscsv="$outdir/${prefix}_B_values.csv"
      if [[ -z "$f16" ]]; then f16="$R16_DIR/merged_files/$cond/__missing_${prefix}.csv"; fi
      if [[ -z "$f24" ]]; then f24="$R24_DIR/merged_files/$cond/__missing_${prefix}.csv"; fi
      echo "[analysis] $cond | $ph | power-law fit"
      run_root "lum_vs_v_systematics.C(\"$f16\",\"$f24\",\"$ph\",\"$syscsv\")"
      run_root "lum_vs_v_fit.C(\"$all20\",\"$ph\",\"$syscsv\",\"$outdir\",\"$prefix\")"
    else
      syscsv="$outdir/${prefix}_lambda_values.csv"
      if [[ -z "$f16" ]]; then f16="$R16_DIR/merged_files/$cond/__missing_${prefix}.csv"; fi
      if [[ -z "$f24" ]]; then f24="$R24_DIR/merged_files/$cond/__missing_${prefix}.csv"; fi
      echo "[analysis] $cond | $ph | exponential fit"
      run_root "lum_vs_T_systematics.C(\"$f16\",\"$f24\",\"$ph\",\"$syscsv\")"
      run_root "lum_vs_T_fit.C(\"$all20\",\"$ph\",\"$syscsv\",\"$outdir\",\"$prefix\")"
    fi
  done
}

run_phases() {
  local c="$1" cond all20 outdir prefix
  cond="$(condition_name "$c")"
  all20="$(all_phases_file "$c")"
  outdir="$MACRO_DIR/phases/$cond"
  prefix="$cond"

  if [[ ! -f "$all20" ]]; then
    echo "WARNING: R=20 all-phases file not found: $all20" >&2
    return 0
  fi

  # Recreate the output for this sensor/constant to prevent stale phase canvases.
  rm -rf "$outdir"
  mkdir -p "$outdir"

  echo "[phases] $cond"
  if [[ "$c" == "T" ]]; then
    run_root "lum_vs_phase_T_const.C(\"$all20\",\"$outdir\",\"$prefix\")"
  else
    run_root "lum_vs_phase_v_const.C(\"$all20\",\"$outdir\",\"$prefix\")"
  fi
}

for cmd in "${COMMANDS[@]}"; do
  for c in "${CONSTANTS[@]}"; do
    if [[ "$cmd" == "analysis" ]]; then run_analysis "$c"; else run_phases "$c"; fi
  done
done

echo "Completed 5analysis.sh for sensor $SENSOR."
