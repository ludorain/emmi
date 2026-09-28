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

usage() {
  cat <<USAGE
Usage: $0 SENSOR [T|v] [analysis|phases]
       $0 SENSOR [analysis|phases]
USAGE
}

if [[ $# -eq 0 ]]; then echo "ERROR: at least the sensor must be specified." >&2; usage >&2; exit 2; fi
SENSOR="$1"; shift
case "$SENSOR" in A1|A2|B1|B2) ;; *) echo "ERROR: invalid sensor '$SENSOR'." >&2; exit 2;; esac

CONSTANT=""; COMMAND=""
for arg in "$@"; do
  case "$arg" in
    T|v) [[ -z "$CONSTANT" ]] || { echo "ERROR: constant specified more than once." >&2; exit 2; }; CONSTANT="$arg" ;;
    analysis|phases) [[ -z "$COMMAND" ]] || { echo "ERROR: command specified more than once." >&2; exit 2; }; COMMAND="$arg" ;;
    *) echo "ERROR: unknown argument '$arg'." >&2; usage >&2; exit 2 ;;
  esac
done
[[ -z "$CONSTANT" ]] && CONSTANTS=(T v) || CONSTANTS=("$CONSTANT")
[[ -z "$COMMAND" ]] && COMMANDS=(analysis phases) || COMMANDS=("$COMMAND")

for f in analysis_common.h phase_common.h lum_vs_v_systematics.C lum_vs_T_systematics.C lum_vs_v_fit.C lum_vs_T_fit.C lum_vs_phase_T_const.C lum_vs_phase_v_const.C; do
  [[ -f "$MACRO_DIR/$f" ]] || { echo "ERROR: missing macro/helper $MACRO_DIR/$f" >&2; exit 3; }
done
command -v root >/dev/null 2>&1 || { echo "ERROR: ROOT executable 'root' was not found in PATH." >&2; exit 3; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 was not found in PATH." >&2; exit 3; }

condition_name() {
  local c="$1"
  if [[ "$c" == "T" ]]; then printf '%s_T=20' "$SENSOR"; else printf '%s_v=5' "$SENSOR"; fi
}

all_phases_file_at_radius() {
  local base="$1" c="$2" cond
  cond="$(condition_name "$c")"
  printf '%s/merged_files/%s/%s_%s_all_phases_global_ID.csv' "$base" "$cond" "$SENSOR" "$c"
}

# Read phases directly from the authoritative R=20 master. This avoids a hard-
# coded phase list silently omitting a valid dataset such as 150 C / 25 h.
list_phases_from_csv() {
  local csv="$1"
  python3 - "$csv" <<'PY'
import csv,re,sys
path=sys.argv[1]
with open(path,newline='',encoding='utf-8-sig') as f:
    r=csv.DictReader(f)
    phases={row.get('phase','').strip() for row in r if row.get('phase','').strip()}
def key(p):
    if p=='before_annealing': return (0,-1e99,-1e99,p)
    m=re.fullmatch(r'annealing_T=([-+]?(?:\d+(?:\.\d*)?|\.\d+))_h=([-+]?(?:\d+(?:\.\d*)?|\.\d+))',p)
    if m: return (1,float(m.group(1)),float(m.group(2)),p)
    return (2,1e99,1e99,p)
for p in sorted(phases,key=key): print(p)
PY
}

run_root() { local expr="$1"; ( cd "$MACRO_DIR" && root -l -b -q "$expr" ); }

require_deltaL_column() {
  local csv="$1"
  python3 - "$csv" <<'PY'
import csv, sys
path=sys.argv[1]
with open(path, newline='', encoding='utf-8-sig') as f:
    r=csv.reader(f)
    header=next(r, [])
if 'deltaL' not in header:
    raise SystemExit(
        f"ERROR: {path} has no deltaL column. Run 4luminosity_sistematical_uncertainty.sh first."
    )
PY
}

run_analysis() {
  local c="$1" cond all20 all16 all24 analysis_base
  cond="$(condition_name "$c")"
  all20="$(all_phases_file_at_radius "$R20_DIR" "$c")"
  all16="$(all_phases_file_at_radius "$R16_DIR" "$c")"
  all24="$(all_phases_file_at_radius "$R24_DIR" "$c")"
  analysis_base="$MACRO_DIR/analysis_chi<4/$cond"

  [[ -f "$all20" ]] || { echo "WARNING: nominal R=20 all-phases file not found: $all20" >&2; return 0; }
  require_deltaL_column "$all20" || return 1
  local phases=() ph
  while IFS= read -r ph; do
    [[ -n "$ph" ]] && phases+=("$ph")
  done < <(list_phases_from_csv "$all20")
  (( ${#phases[@]} > 0 )) || { echo "WARNING: no phases found in $all20" >&2; return 0; }

  for ph in "${phases[@]}"; do
    local outdir prefix syscsv
    outdir="$analysis_base/$ph"; prefix="${cond}_${ph}"
    rm -rf "$outdir"; mkdir -p "$outdir"

    if [[ "$c" == "T" ]]; then
      syscsv="$outdir/${prefix}_B_values.csv"
      echo "[analysis] $cond | $ph | power-law fit"
      if [[ -f "$all16" && -f "$all24" ]]; then
        run_root "lum_vs_v_systematics.C(\"$all20\",\"$all16\",\"$all24\",\"$ph\",\"$syscsv\")"
      else
        echo "WARNING: R16/R24 master missing: parameter systematic for B will be unavailable." >&2
        : > "$syscsv"
      fi
      run_root "lum_vs_v_fit.C(\"$all20\",\"$ph\",\"$syscsv\",\"$outdir\",\"$prefix\")"
    else
      syscsv="$outdir/${prefix}_lambda_values.csv"
      echo "[analysis] $cond | $ph | exponential fit"
      if [[ -f "$all16" && -f "$all24" ]]; then
        run_root "lum_vs_T_systematics.C(\"$all20\",\"$all16\",\"$all24\",\"$ph\",\"$syscsv\")"
      else
        echo "WARNING: R16/R24 master missing: parameter systematic for lambda will be unavailable." >&2
        : > "$syscsv"
      fi
      run_root "lum_vs_T_fit.C(\"$all20\",\"$ph\",\"$syscsv\",\"$outdir\",\"$prefix\")"
    fi
  done
}

run_phases() {
  local c="$1" cond all20 all16 all24 outdir prefix
  cond="$(condition_name "$c")"
  all20="$(all_phases_file_at_radius "$R20_DIR" "$c")"
  all16="$(all_phases_file_at_radius "$R16_DIR" "$c")"
  all24="$(all_phases_file_at_radius "$R24_DIR" "$c")"
  outdir="$MACRO_DIR/phases/$cond"; prefix="$cond"
  [[ -f "$all20" ]] || { echo "WARNING: R=20 all-phases file not found: $all20" >&2; return 0; }
  require_deltaL_column "$all20" || return 1
  [[ -f "$all16" ]] || all16=""
  [[ -f "$all24" ]] || all24=""
  rm -rf "$outdir"; mkdir -p "$outdir"

  echo "[phases] $cond"
  if [[ "$c" == "T" ]]; then
    run_root "lum_vs_phase_T_const.C(\"$all20\",\"$outdir\",\"$prefix\",\"$all16\",\"$all24\")"
  else
    run_root "lum_vs_phase_v_const.C(\"$all20\",\"$outdir\",\"$prefix\",\"$all16\",\"$all24\")"
  fi
}

for cmd in "${COMMANDS[@]}"; do
  for c in "${CONSTANTS[@]}"; do
    if [[ "$cmd" == "analysis" ]]; then run_analysis "$c"; else run_phases "$c"; fi
  done
done

echo "Completed 5analysis.sh for sensor $SENSOR."
