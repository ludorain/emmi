#!/usr/bin/env bash

# =============================================================================
# 2analysis.sh
# Reduced-sample EMMI analysis for emmi/new_device/gold
#
# Input produced by 1file_processing.sh:
#   DATA/merged_files/A1_all_runs_global_ID.csv          (nominal R=20 + deltaL)
#   DATA_R=16/merged_files/A1_all_runs_global_ID.csv
#   DATA_R=24/merged_files/A1_all_runs_global_ID.csv
#
# Analyses:
#   T20 : L(Vover)=A*Vover^B, per-hotspot luminosity systematics, deltaB,
#         B vs global spot ID with statistical+systematic errors and B=2 line.
#   v7  : L(T)=A*exp(lambda*T), per-hotspot luminosity systematics, deltaLambda,
#         lambda vs global spot ID + stat-only constant fit.
#   v9  : same as v7.
#   v7/v9 comparison: overlaid L vs T per hotspot and overlaid lambda vs spot.
#
# No annealing-phase option is present in this reduced pipeline.
# =============================================================================

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO_DIR="${MACRO_DIR:-$SCRIPT_DIR/analysis_macros}"
OUTPUT_DIR="${OUTPUT_DIR:-$SCRIPT_DIR/analysis}"

R20_CSV="${R20_CSV:-$SCRIPT_DIR/DATA/merged_files/A1_all_runs_global_ID.csv}"
R16_CSV="${R16_CSV:-$SCRIPT_DIR/DATA_R=16/merged_files/A1_all_runs_global_ID.csv}"
R24_CSV="${R24_CSV:-$SCRIPT_DIR/DATA_R=24/merged_files/A1_all_runs_global_ID.csv}"

fail(){ echo "ERROR: $*" >&2; exit 1; }
run_root(){ local expr="$1"; ( cd "$MACRO_DIR" && root -l -b -q "$expr" ); }

check_csv(){
    local f="$1"
    [[ -f "$f" ]] || fail "required CSV not found: $f"
}

require_columns(){
    local f="$1"; shift
    python3 - "$f" "$@" <<'PY'
import csv,sys
path=sys.argv[1]; required=sys.argv[2:]
with open(path,newline='',encoding='utf-8-sig') as h:
    r=csv.reader(h); header=next(r,[])
missing=[c for c in required if c not in header]
if missing:
    raise SystemExit(f"ERROR: {path} is missing columns: {missing}")
PY
}

main(){
    command -v root >/dev/null 2>&1 || fail "ROOT executable 'root' not found in PATH"
    command -v python3 >/dev/null 2>&1 || fail "python3 not found in PATH"

    for f in gold_analysis_common.h gold_v_systematics.C gold_T_systematics.C gold_lum_vs_v_fit.C gold_lum_vs_T_fit.C gold_compare_v7_v9.C; do
        [[ -f "$MACRO_DIR/$f" ]] || fail "missing macro/helper: $MACRO_DIR/$f"
    done

    check_csv "$R20_CSV"; check_csv "$R16_CSV"; check_csv "$R24_CSV"
    require_columns "$R20_CSV" spot x y luminosity error T v v_fin dataset_key detected deltaL
    require_columns "$R16_CSV" spot x y luminosity error T v v_fin dataset_key detected
    require_columns "$R24_CSV" spot x y luminosity error T v v_fin dataset_key detected

    rm -rf "$OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR/A1_T=20" "$OUTPUT_DIR/A1_v=7" "$OUTPUT_DIR/A1_v=9" "$OUTPUT_DIR/comparison_v7_v9"

    echo "============================================================"
    echo "GOLD SAMPLE ANALYSIS"
    echo "R=20:  $R20_CSV"
    echo "R=16:  $R16_CSV"
    echo "R=24:  $R24_CSV"
    echo "Output: $OUTPUT_DIR"
    echo "============================================================"

    # -------------------------------------------------------------------------
    # A1_T=20: power-law analysis
    # -------------------------------------------------------------------------
    local B_SYS="$OUTPUT_DIR/A1_T=20/A1_T=20_B_systematics.csv"
    echo "[1/6] A1_T=20 - systematic uncertainty on B"
    run_root "gold_v_systematics.C(\"$R20_CSV\",\"$R16_CSV\",\"$R24_CSV\",\"$B_SYS\")"

    echo "[2/6] A1_T=20 - luminosity vs overvoltage + B vs spot"
    run_root "gold_lum_vs_v_fit.C(\"$R20_CSV\",\"$B_SYS\",\"$OUTPUT_DIR/A1_T=20\")"

    # -------------------------------------------------------------------------
    # A1_v=7 and A1_v=9: exponential temperature analyses
    # -------------------------------------------------------------------------
    local L7_SYS="$OUTPUT_DIR/A1_v=7/A1_v=7_lambda_systematics.csv"
    local L9_SYS="$OUTPUT_DIR/A1_v=9/A1_v=9_lambda_systematics.csv"

    echo "[3/6] A1_v=7 - lambda systematics + luminosity vs T"
    run_root "gold_T_systematics.C(\"$R20_CSV\",\"$R16_CSV\",\"$R24_CSV\",\"v7\",\"$L7_SYS\")"
    run_root "gold_lum_vs_T_fit.C(\"$R20_CSV\",\"v7\",\"$L7_SYS\",\"$OUTPUT_DIR/A1_v=7\")"

    echo "[4/6] A1_v=9 - lambda systematics + luminosity vs T"
    run_root "gold_T_systematics.C(\"$R20_CSV\",\"$R16_CSV\",\"$R24_CSV\",\"v9\",\"$L9_SYS\")"
    run_root "gold_lum_vs_T_fit.C(\"$R20_CSV\",\"v9\",\"$L9_SYS\",\"$OUTPUT_DIR/A1_v=9\")"

    # -------------------------------------------------------------------------
    # Combined v=7 / v=9 canvases
    # -------------------------------------------------------------------------
    local L7_RESULTS="$OUTPUT_DIR/A1_v=7/A1_v7_lambda_fit_results.csv"
    local L9_RESULTS="$OUTPUT_DIR/A1_v=9/A1_v9_lambda_fit_results.csv"

    [[ -f "$L7_RESULTS" ]] || fail "lambda result file not produced: $L7_RESULTS"
    [[ -f "$L9_RESULTS" ]] || fail "lambda result file not produced: $L9_RESULTS"

    echo "[5/6] v=7 / v=9 - overlaid luminosity vs temperature"
    echo "[6/6] v=7 / v=9 - overlaid lambda vs spot"
    run_root "gold_compare_v7_v9.C(\"$R20_CSV\",\"$L7_RESULTS\",\"$L9_RESULTS\",\"$OUTPUT_DIR/comparison_v7_v9\")"

    echo
    echo "============================================================"
    echo "2analysis.sh COMPLETED SUCCESSFULLY"
    echo "Output directory: $OUTPUT_DIR"
    echo "============================================================"
}

main "$@"
