#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 6luminosity_distribution.sh
#
# Build luminosity distributions at T = 20 C, using the requested
# overvoltage point for each irradiated sensor:
#     A1: v_fin = 7 V
#     B1: v_fin = 5 V
# for:
#   - irradiated A1/B1 sensors, for every annealing phase;
#   - the A1 new-device reference sample.
#
# The script first extracts the requested rows/columns into small,
# auditable CSV files, then passes those files to the ROOT macro
# luminosity_distribution.C.
#
# Place this script and luminosity_distribution.C in emmi/ and run:
#     ./6luminosity_distribution.sh
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ANALYSIS_DIR="$SCRIPT_DIR/spot_luminosity_final/analysis"
NEW_DEVICE_CSV="$SCRIPT_DIR/new_device/gold/DATA/A1_v=7_run=20250826-035836/6luminosity/A1_v=7_run=20250826-035836_global_ID.csv"
MACRO="$SCRIPT_DIR/luminosity_distribution.C"

OUTPUT_DIR="$SCRIPT_DIR/luminosity_distributions"
SELECTED_DIR="$OUTPUT_DIR/selected_csv"

SENSORS=(A1 B1)
PHASES=(
    before_annealing
    annealing_T=75_h=5
    annealing_T=75_h=25
    annealing_T=100_h=5
    annealing_T=100_h=25
    annealing_T=125_h=5
    annealing_T=125_h=25
    annealing_T=150_h=5
    annealing_T=150_h=25
)

mkdir -p "$OUTPUT_DIR" "$SELECTED_DIR"

if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 not found." >&2
    exit 1
fi

if ! command -v root >/dev/null 2>&1; then
    echo "ERROR: ROOT executable 'root' not found." >&2
    exit 1
fi

if [[ ! -f "$MACRO" ]]; then
    echo "ERROR: ROOT macro not found: $MACRO" >&2
    exit 1
fi

# ------------------------------------------------------------
# Extract one selected sample.
#
# kind = irradiated:
#   keep v_fin = target_vfin and detection = present
#
# kind = new:
#   keep T = 20, v = 7 and detected = True
#
# Output columns:
#   luminosity,error,deltaL
# ------------------------------------------------------------
extract_selected_csv() {
    local input_csv="$1"
    local output_csv="$2"
    local kind="$3"
    local target_vfin="${4:-}"

    if [[ ! -f "$input_csv" ]]; then
        echo "ERROR: input CSV not found: $input_csv" >&2
        exit 1
    fi

    python3 - "$input_csv" "$output_csv" "$kind" "$target_vfin" <<'PY'
import csv
import math
import os
import sys

input_csv, output_csv, kind, target_vfin_arg = sys.argv[1:5]


def detect_delimiter(path):
    with open(path, "r", newline="", encoding="utf-8-sig") as f:
        first = f.readline()
    return ";" if first.count(";") > first.count(",") else ","


def clean_header(name):
    return name.strip().lstrip("\ufeff")


def parse_float(value, column, line_number):
    try:
        x = float(str(value).strip())
    except Exception as exc:
        raise RuntimeError(
            f"invalid numeric value in column '{column}' at line {line_number}: {value!r}"
        ) from exc
    if not math.isfinite(x):
        raise RuntimeError(
            f"non-finite value in column '{column}' at line {line_number}: {value!r}"
        )
    return x


delimiter = detect_delimiter(input_csv)

with open(input_csv, "r", newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter=delimiter)
    if reader.fieldnames is None:
        raise RuntimeError(f"CSV has no header: {input_csv}")

    # Re-map possible whitespace/BOM in the original header.
    field_map = {clean_header(name): name for name in reader.fieldnames}

    required = {"luminosity", "error", "deltaL"}
    if kind == "irradiated":
        required |= {"v_fin", "detection"}
        if not target_vfin_arg:
            raise RuntimeError("missing target v_fin for irradiated selection")
        target_vfin = float(target_vfin_arg)
    elif kind == "new":
        required |= {"T", "v", "detected"}
    else:
        raise RuntimeError(f"unknown extraction kind: {kind}")

    missing = sorted(required - set(field_map))
    if missing:
        raise RuntimeError(
            f"missing column(s) {missing} in {input_csv}. "
            f"Available columns: {sorted(field_map)}"
        )

    selected = []
    rejected_bad = 0

    for line_number, row in enumerate(reader, start=2):
        try:
            if kind == "irradiated":
                v_fin = parse_float(row[field_map["v_fin"]], "v_fin", line_number)
                detection = str(row[field_map["detection"]]).strip().lower()

                if abs(v_fin - target_vfin) > 1.0e-6:
                    continue
                if detection != "present":
                    continue

            else:  # new device
                temperature = parse_float(row[field_map["T"]], "T", line_number)
                voltage = parse_float(row[field_map["v"]], "v", line_number)
                detected = str(row[field_map["detected"]]).strip().lower()

                if abs(temperature - 20.0) > 1.0e-6:
                    continue
                if abs(voltage - 7.0) > 1.0e-6:
                    continue
                if detected not in {"true", "1", "yes"}:
                    continue

            luminosity = parse_float(
                row[field_map["luminosity"]], "luminosity", line_number
            )
            stat_error = parse_float(
                row[field_map["error"]], "error", line_number
            )
            syst_error = parse_float(
                row[field_map["deltaL"]], "deltaL", line_number
            )

            selected.append((luminosity, stat_error, syst_error))

        except RuntimeError as exc:
            rejected_bad += 1
            print(f"WARNING: {exc}", file=sys.stderr)

os.makedirs(os.path.dirname(output_csv), exist_ok=True)
with open(output_csv, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["luminosity", "error", "deltaL"])
    writer.writerows(selected)

print(
    f"Selected {len(selected)} rows from {input_csv} -> {output_csv}"
    + (f" ({rejected_bad} malformed row(s) skipped)" if rejected_bad else "")
)

if not selected:
    raise RuntimeError(
        f"selection produced zero rows for {input_csv} (kind={kind})"
    )
PY
}

# ROOT command helper. The paths used in this project do not contain
# double quotes; escape backslashes/double quotes anyway for safety.
cxx_escape() {
    local s="$1"
    s="${s//\\/\\\\}"
    s="${s//\"/\\\"}"
    printf '%s' "$s"
}

run_root() {
    local action="$1"
    local csv1="$2"
    local csv2="$3"
    local sensor="$4"
    local phase="$5"

    local macro_e csv1_e csv2_e out_e sensor_e phase_e action_e
    macro_e="$(cxx_escape "$MACRO")"
    csv1_e="$(cxx_escape "$csv1")"
    csv2_e="$(cxx_escape "$csv2")"
    out_e="$(cxx_escape "$OUTPUT_DIR")"
    sensor_e="$(cxx_escape "$sensor")"
    phase_e="$(cxx_escape "$phase")"
    action_e="$(cxx_escape "$action")"

    root -l -b -q \
        "${macro_e}(\"${action_e}\",\"${csv1_e}\",\"${csv2_e}\",\"${sensor_e}\",\"${phase_e}\",\"${out_e}\")"
}

# ============================================================
# 1. Extract irradiated samples and produce individual plots.
# ============================================================

echo
printf '%s\n' "============================================================"
printf '%s\n' "Selecting irradiated samples: A1 v_fin = 7, B1 v_fin = 5, detection = present"
printf '%s\n' "============================================================"

for sensor in "${SENSORS[@]}"; do
    case "$sensor" in
        A1) target_vfin="7" ;;
        B1) target_vfin="5" ;;
        *)
            echo "ERROR: unsupported sensor for v_fin selection: $sensor" >&2
            exit 1
            ;;
    esac

    for phase in "${PHASES[@]}"; do
        source_csv="$ANALYSIS_DIR/${sensor}_T=20/${phase}/${sensor}_T=20_${phase}_analysis.csv"
        selected_csv="$SELECTED_DIR/${sensor}_${phase}_T20_vfin${target_vfin}_present.csv"

        extract_selected_csv "$source_csv" "$selected_csv" irradiated "$target_vfin"
        run_root individual "$selected_csv" "" "$sensor" "$phase"
    done
done

# ============================================================
# 2. Extract A1 new-device reference at T = 20 C.
# ============================================================

echo
printf '%s\n' "============================================================"
printf '%s\n' "Selecting A1 new-device sample: T = 20 C, v = 7, detected = True"
printf '%s\n' "============================================================"

NEW_SELECTED="$SELECTED_DIR/A1_new_device_T20_v7.csv"
extract_selected_csv "$NEW_DEVICE_CSV" "$NEW_SELECTED" new

# ============================================================
# Produce also the non-normalized luminosity distribution for A1 new device.
# ============================================================
run_root individual "$NEW_SELECTED" "" A1 new_device

# ============================================================
# 3. A1 new device vs A1 before annealing.
#    The two distributions are normalized independently and overlaid.
#    No bin-by-bin ratio is produced.
# ============================================================

A1_BEFORE_SELECTED="$SELECTED_DIR/A1_before_annealing_T20_vfin7_present.csv"
run_root new_before \
    "$NEW_SELECTED" \
    "$A1_BEFORE_SELECTED" \
    A1 \
    before_annealing

# ============================================================
# 4. Before annealing vs each subsequent annealing phase.
#    The two distributions are normalized independently and overlaid.
#    No bin-by-bin ratio is produced.
# ============================================================

for sensor in "${SENSORS[@]}"; do
    case "$sensor" in
        A1) target_vfin="7" ;;
        B1) target_vfin="5" ;;
        *)
            echo "ERROR: unsupported sensor for v_fin selection: $sensor" >&2
            exit 1
            ;;
    esac

    before_selected="$SELECTED_DIR/${sensor}_before_annealing_T20_vfin${target_vfin}_present.csv"

    for phase in "${PHASES[@]:1}"; do
        after_selected="$SELECTED_DIR/${sensor}_${phase}_T20_vfin${target_vfin}_present.csv"

        run_root phase_compare \
            "$before_selected" \
            "$after_selected" \
            "$sensor" \
            "$phase"
    done
done

echo
echo "Done."
echo "Plots saved in: $OUTPUT_DIR"
echo "Selected input samples saved in: $SELECTED_DIR"
