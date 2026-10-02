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
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MACRO_DIR="$SCRIPT_DIR/spot_luminosity_final"

R16_DIR="$SCRIPT_DIR/DATA_irradiated_isolated_R=16"
R20_DIR="$SCRIPT_DIR/DATA_irradiated_isolated_R=20"
R24_DIR="$SCRIPT_DIR/DATA_irradiated_isolated_R=24"


# ============================================================
# Usage
# ============================================================

usage() {
  cat <<USAGE
Usage: $0 SENSOR [T|v] [analysis|phases]
       $0 SENSOR [analysis|phases]
USAGE
}


# ============================================================
# Parse sensor
# ============================================================

if [[ $# -eq 0 ]]; then
  echo "ERROR: at least the sensor must be specified." >&2
  usage >&2
  exit 2
fi

SENSOR="$1"
shift

case "$SENSOR" in
  A1|A2|B1|B2)
    ;;
  *)
    echo "ERROR: invalid sensor '$SENSOR'." >&2
    exit 2
    ;;
esac


# ============================================================
# Parse optional arguments
# ============================================================

CONSTANT=""
COMMAND=""

for arg in "$@"; do

  case "$arg" in

    T|v)
      [[ -z "$CONSTANT" ]] || {
        echo "ERROR: constant specified more than once." >&2
        exit 2
      }
      CONSTANT="$arg"
      ;;

    analysis|phases)
      [[ -z "$COMMAND" ]] || {
        echo "ERROR: command specified more than once." >&2
        exit 2
      }
      COMMAND="$arg"
      ;;

    *)
      echo "ERROR: unknown argument '$arg'." >&2
      usage >&2
      exit 2
      ;;

  esac

done


# If no constant is specified -> run both T and v.
if [[ -z "$CONSTANT" ]]; then
  CONSTANTS=(T v)
else
  CONSTANTS=("$CONSTANT")
fi


# If no command is specified -> run both analysis and phases.
if [[ -z "$COMMAND" ]]; then
  COMMANDS=(analysis phases)
else
  COMMANDS=("$COMMAND")
fi


# ============================================================
# Check required files/programs
# ============================================================

for f in \
  analysis_common.h \
  phase_common.h \
  lum_vs_v_systematics.C \
  lum_vs_T_systematics.C \
  lum_vs_v_fit.C \
  lum_vs_T_fit.C \
  lum_vs_phase_T_const.C \
  lum_vs_phase_v_const.C
do

  [[ -f "$MACRO_DIR/$f" ]] || {
    echo "ERROR: missing macro/helper $MACRO_DIR/$f" >&2
    exit 3
  }

done


command -v root >/dev/null 2>&1 || {
  echo "ERROR: ROOT executable 'root' was not found in PATH." >&2
  exit 3
}

command -v python3 >/dev/null 2>&1 || {
  echo "ERROR: python3 was not found in PATH." >&2
  exit 3
}


# ============================================================
# Condition name
#
# T means:
#   temperature is constant at T = 20 C
#   -> luminosity vs overvoltage
#   -> fit parameter B
#
# v means:
#   voltage/overvoltage is constant
#   -> luminosity vs temperature
#   -> fit parameter lambda
# ============================================================

condition_name() {

  local c="$1"

  if [[ "$c" == "T" ]]; then

    printf '%s_T=20' "$SENSOR"

  else

    case "$SENSOR" in

      B1)
        printf '%s_v=3' "$SENSOR"
        ;;

      *)
        printf '%s_v=5' "$SENSOR"
        ;;

    esac

  fi
}


# ============================================================
# Path of the all-phases master CSV for a given radius
# ============================================================

all_phases_file_at_radius() {

  local base="$1"
  local c="$2"
  local cond

  cond="$(condition_name "$c")"

  printf \
    '%s/merged_files/%s/%s_%s_all_phases_global_ID.csv' \
    "$base" \
    "$cond" \
    "$SENSOR" \
    "$c"
}


# ============================================================
# Read annealing phases from the authoritative R=20 master CSV
# ============================================================

list_phases_from_csv() {

  local csv="$1"

  python3 - "$csv" <<'PY'
import csv
import re
import sys

path = sys.argv[1]

with open(path, newline='', encoding='utf-8-sig') as f:
    r = csv.DictReader(f)

    phases = {
        row.get('phase', '').strip()
        for row in r
        if row.get('phase', '').strip()
    }

def key(p):

    if p == 'before_annealing':
        return (0, -1e99, -1e99, p)

    m = re.fullmatch(
        r'annealing_T=([-+]?(?:\d+(?:\.\d*)?|\.\d+))_h=([-+]?(?:\d+(?:\.\d*)?|\.\d+))',
        p
    )

    if m:
        return (
            1,
            float(m.group(1)),
            float(m.group(2)),
            p
        )

    return (2, 1e99, 1e99, p)


for p in sorted(phases, key=key):
    print(p)

PY
}


# ============================================================
# ROOT runner
# ============================================================

run_root() {

  local expr="$1"

  (
    cd "$MACRO_DIR" &&
    root -l -b -q "$expr"
  )
}


# ============================================================
# Check that the R=20 master contains deltaL
# ============================================================

require_deltaL_column() {

  local csv="$1"

  python3 - "$csv" <<'PY'
import csv
import sys

path = sys.argv[1]

with open(path, newline='', encoding='utf-8-sig') as f:

    r = csv.reader(f)
    header = next(r, [])

if 'deltaL' not in header:

    raise SystemExit(
        f"ERROR: {path} has no deltaL column. "
        "Run 4luminosity_sistematical_uncertainty.sh first."
    )

PY
}


# ============================================================
# ANALYSIS
# ============================================================

run_analysis() {

  local c="$1"

  local cond
  local all20
  local all16
  local all24
  local analysis_base

  # Summary files:
  #
  # lambda_summary_csv:
  #   phase, lambda_mean, lambda_mean_stat_error, lambda_rms, n_spots
  #
  # B_summary_csv:
  #   reserved for the analogous B-vs-phase summary.
  #
  # No "global" terminology is used anymore.

  local lambda_summary_csv
  local B_summary_csv


  # ----------------------------------------------------------
  # Resolve paths
  # ----------------------------------------------------------

  cond="$(condition_name "$c")"

  all20="$(all_phases_file_at_radius "$R20_DIR" "$c")"
  all16="$(all_phases_file_at_radius "$R16_DIR" "$c")"
  all24="$(all_phases_file_at_radius "$R24_DIR" "$c")"

  analysis_base="$MACRO_DIR/analysis/$cond"

  mkdir -p "$analysis_base"


  # ----------------------------------------------------------
  # Prepare phase-summary CSV
  # ----------------------------------------------------------

  lambda_summary_csv=""
  B_summary_csv=""


  if [[ "$c" == "T" ]]; then

    # --------------------------------------------------------
    # T constant -> Luminosity vs V -> parameter B
    #
    # The B code will later be modified analogously to lambda.
    # We already use the neutral/non-global filename here.
    # --------------------------------------------------------

    B_summary_csv="$analysis_base/${cond}_B_vs_phase.csv"

    # Start from a clean summary file for this analysis run.
    rm -f "$B_summary_csv"

  else

    # --------------------------------------------------------
    # v constant -> Luminosity vs T -> parameter lambda
    #
    # The file will contain, for every annealing phase:
    #
    # phase
    # lambda_mean
    # lambda_mean_stat_error
    # lambda_rms
    # n_spots
    # --------------------------------------------------------

    lambda_summary_csv="$analysis_base/${cond}_lambda_vs_phase.csv"

    # Start from a clean summary file for this analysis run.
    # lum_vs_T_fit.C will append one row per phase.
    rm -f "$lambda_summary_csv"

  fi


  # ----------------------------------------------------------
  # Check nominal master
  # ----------------------------------------------------------

  [[ -f "$all20" ]] || {

    echo \
      "WARNING: nominal R=20 all-phases file not found: $all20" \
      >&2

    return 0
  }


  require_deltaL_column "$all20" || return 1


  # ----------------------------------------------------------
  # Read available annealing phases
  # ----------------------------------------------------------

  local phases=()
  local ph

  while IFS= read -r ph; do

    [[ -n "$ph" ]] &&
      phases+=("$ph")

  done < <(list_phases_from_csv "$all20")


  (( ${#phases[@]} > 0 )) || {

    echo \
      "WARNING: no phases found in $all20" \
      >&2

    return 0
  }


  # ==========================================================
  # Loop over annealing phases
  # ==========================================================

  for ph in "${phases[@]}"; do

    local outdir
    local prefix
    local syscsv

    outdir="$analysis_base/$ph"
    prefix="${cond}_${ph}"

    rm -rf "$outdir"
    mkdir -p "$outdir"


    # ========================================================
    # T constant
    #
    # -> Luminosity vs overvoltage
    # -> power-law exponent B
    # ========================================================

    if [[ "$c" == "T" ]]; then

      syscsv="$outdir/${prefix}_B_values.csv"

      echo \
        "[analysis] $cond | $ph | power-law fit"


      # ------------------------------------------------------
      # Systematic uncertainty on B
      # ------------------------------------------------------

      if [[ -f "$all16" && -f "$all24" ]]; then

        run_root \
          "lum_vs_v_systematics.C(\"$all20\",\"$all16\",\"$all24\",\"$ph\",\"$syscsv\")"

      else

        echo \
          "WARNING: R16/R24 master missing: parameter systematic for B will be unavailable." \
          >&2

        : > "$syscsv"

      fi


      # ------------------------------------------------------
      # Nominal B analysis
      #
      # The last argument is the B-vs-phase summary CSV.
      #
      # When lum_vs_v_fit.C is updated analogously to lambda,
      # this file will contain mean B, statistical error on
      # the mean, RMS and number of hotspots.
      # ------------------------------------------------------

      run_root \
        "lum_vs_v_fit.C(\"$all20\",\"$ph\",\"$syscsv\",\"$outdir\",\"$prefix\",\"\",\"\",\"$B_summary_csv\")"


    # ========================================================
    # v constant
    #
    # -> Luminosity vs temperature
    # -> exponential parameter lambda
    # ========================================================

    else

      syscsv="$outdir/${prefix}_lambda_values.csv"

      echo \
        "[analysis] $cond | $ph | exponential fit"


      # ------------------------------------------------------
      # Systematic uncertainty on lambda for each hotspot
      # ------------------------------------------------------

      if [[ -f "$all16" && -f "$all24" ]]; then

        run_root \
          "lum_vs_T_systematics.C(\"$all20\",\"$all16\",\"$all24\",\"$ph\",\"$syscsv\")"

      else

        echo \
          "WARNING: R16/R24 master missing: parameter systematic for lambda will be unavailable." \
          >&2

        : > "$syscsv"

      fi


      # ------------------------------------------------------
      # Nominal lambda analysis
      #
      # lum_vs_T_fit.C now:
      #
      # 1. fits every hotspot independently;
      # 2. selects acceptable lambda_i values;
      # 3. calculates:
      #       lambda_mean
      #       lambda_mean_stat_error
      #       lambda_rms
      # 4. appends one row to lambda_summary_csv.
      #
      # The two empty arguments correspond to the optional
      # R16/R24 comparison files, which are not passed here.
      # ------------------------------------------------------

      run_root \
        "lum_vs_T_fit.C(\"$all20\",\"$ph\",\"$syscsv\",\"$outdir\",\"$prefix\",\"\",\"\",\"$lambda_summary_csv\")"

    fi

  done
}


# ============================================================
# PHASE EVOLUTION
# ============================================================

run_phases() {

  local c="$1"

  local cond
  local all20
  local all16
  local all24
  local outdir
  local prefix


  cond="$(condition_name "$c")"

  all20="$(all_phases_file_at_radius "$R20_DIR" "$c")"
  all16="$(all_phases_file_at_radius "$R16_DIR" "$c")"
  all24="$(all_phases_file_at_radius "$R24_DIR" "$c")"

  outdir="$MACRO_DIR/phases/$cond"
  prefix="$cond"


  [[ -f "$all20" ]] || {

    echo \
      "WARNING: R=20 all-phases file not found: $all20" \
      >&2

    return 0
  }


  require_deltaL_column "$all20" || return 1


  # R16/R24 are optional in phase-comparison mode.
  [[ -f "$all16" ]] || all16=""
  [[ -f "$all24" ]] || all24=""


  rm -rf "$outdir"
  mkdir -p "$outdir"


  echo "[phases] $cond"


  # ----------------------------------------------------------
  # T constant -> luminosity vs phase at fixed T
  # ----------------------------------------------------------

  if [[ "$c" == "T" ]]; then

    run_root \
      "lum_vs_phase_T_const.C(\"$all20\",\"$outdir\",\"$prefix\",\"$all16\",\"$all24\")"


  # ----------------------------------------------------------
  # v constant -> luminosity vs phase at fixed v
  # ----------------------------------------------------------

  else

    run_root \
      "lum_vs_phase_v_const.C(\"$all20\",\"$outdir\",\"$prefix\",\"$all16\",\"$all24\")"

  fi
}


# ============================================================
# Main execution
# ============================================================

for cmd in "${COMMANDS[@]}"; do

  for c in "${CONSTANTS[@]}"; do

    if [[ "$cmd" == "analysis" ]]; then

      run_analysis "$c"

    else

      run_phases "$c"

    fi

  done

done


echo "Completed 5analysis.sh for sensor $SENSOR."