#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# hotspot_phase_gold_error.sh
#
# Read GOLD-region phase-change summaries for A1 and B1
# at pixel sizes 24, 30 and 36.
#
# Nominal dataset:
#   pixel = 30
#
# Systematic uncertainty for each quantity:
#
#   deltaX1 = |X36 - X30|
#   deltaX2 = |X30 - X24|
#   deltaX  = max(deltaX1, deltaX2)
#
# Outputs:
#
#   hotspot_phase_summary_gold_error.csv
#
#       sensor
#       phase
#       counted_hotspots
#       counted_hotspots_error
#
#
#   hotspot_phase_complete_gold_error.csv
#
#       sensor
#       phase
#       run_number
#       counted_hotspots
#       counted_hotspots_error
#       reference
#       reference_error
#       appeared
#       appeared_error
#       reappeared
#       reappeared_error
#       disappeared
#       disappeared_error
#       absent
#       absent_error
#       not_yet_detected
#       not_yet_detected_error
#
# ============================================================


# ------------------------------------------------------------
# Directory containing this script
# ------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# ------------------------------------------------------------
# emmi directory
# ------------------------------------------------------------
EMMI_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"


# ------------------------------------------------------------
# Input directories
# ------------------------------------------------------------
DIR_24="${EMMI_DIR}/DATA_irradiated_NOisolated_changeR_pixel_24/merged_files"
DIR_30="${EMMI_DIR}/DATA_irradiated_NOisolated_changeR_pixel_30/merged_files"
DIR_36="${EMMI_DIR}/DATA_irradiated_NOisolated_changeR_pixel_36/merged_files"


# ------------------------------------------------------------
# A1 input files
# ------------------------------------------------------------
A1_24="${DIR_24}/A1_T=20/A1_T_phase_changes.csv"
A1_30="${DIR_30}/A1_T=20/A1_T_phase_changes.csv"
A1_36="${DIR_36}/A1_T=20/A1_T_phase_changes.csv"


# ------------------------------------------------------------
# B1 input files
# ------------------------------------------------------------
B1_24="${DIR_24}/B1_T=20/B1_T_phase_changes.csv"
B1_30="${DIR_30}/B1_T=20/B1_T_phase_changes.csv"
B1_36="${DIR_36}/B1_T=20/B1_T_phase_changes.csv"


# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------
OUTPUT_SUMMARY="${SCRIPT_DIR}/hotspot_phase_summary_gold_error.csv"

OUTPUT_COMPLETE="${SCRIPT_DIR}/hotspot_phase_complete_gold_error.csv"


# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------
INPUT_FILES=(
    "$A1_24"
    "$A1_30"
    "$A1_36"
    "$B1_24"
    "$B1_30"
    "$B1_36"
)

for file in "${INPUT_FILES[@]}"; do

    if [[ ! -f "$file" ]]; then
        echo
        echo "ERROR: input file not found:"
        echo "  $file"
        echo
        exit 1
    fi

done


echo
echo "============================================================"
echo " GOLD hotspot phase summary + systematic errors"
echo "============================================================"
echo
echo "Nominal dataset:"
echo "  pixel = 30"
echo
echo "Systematic variations:"
echo "  pixel = 24"
echo "  pixel = 36"
echo
echo "Sensors:"
echo "  A1"
echo "  B1"
echo


python3 - \
    "$A1_24" "$A1_30" "$A1_36" \
    "$B1_24" "$B1_30" "$B1_36" \
    "$OUTPUT_SUMMARY" "$OUTPUT_COMPLETE" <<'PY'

import csv
import sys


# ============================================================
# Arguments
# ============================================================

A1_24 = sys.argv[1]
A1_30 = sys.argv[2]
A1_36 = sys.argv[3]

B1_24 = sys.argv[4]
B1_30 = sys.argv[5]
B1_36 = sys.argv[6]

output_summary = sys.argv[7]
output_complete = sys.argv[8]


# ============================================================
# Input -> output column mapping
# ============================================================

COLUMN_MAP = {
    "n_detected": "counted_hotspots",
    "n_reference": "reference",
    "n_appeared": "appeared",
    "n_reappeared": "reappeared",
    "n_disappeared": "disappeared",
    "n_absent": "absent",
    "n_not_yet_detected": "not_yet_detected",
}


# ============================================================
# Required input columns
# ============================================================

REQUIRED_COLUMNS = {
    "phase",
    "run_number",
    "n_detected",
    "n_reference",
    "n_appeared",
    "n_reappeared",
    "n_disappeared",
    "n_absent",
    "n_not_yet_detected",
}


# ============================================================
# Utility
# ============================================================

def clean(value):

    if value is None:
        return ""

    return str(value).strip()


def parse_integer(value, filename, line_number, column):

    value = clean(value)

    try:
        return int(value)

    except ValueError:

        raise RuntimeError(
            f"Invalid integer value '{value}' "
            f"for column '{column}' "
            f"in {filename}, line {line_number}"
        )


# ============================================================
# Read one phase_changes file
#
# The source file has structure:
#
# phase
# run_number
# n_detected
# n_reference
# n_appeared
# n_reappeared
# n_disappeared
# n_absent
# n_not_yet_detected
#
# Internally we rename them to cleaner output names.
# ============================================================

def read_phase_file(filename, sensor):

    rows = []
    data = {}

    with open(
        filename,
        "r",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        reader = csv.DictReader(f)

        if reader.fieldnames is None:
            raise RuntimeError(
                f"CSV has no header: {filename}"
            )

        reader.fieldnames = [
            name.strip()
            for name in reader.fieldnames
            if name is not None
        ]


        missing = REQUIRED_COLUMNS - set(reader.fieldnames)

        if missing:

            raise RuntimeError(
                f"Missing required columns in {filename}: "
                + ", ".join(sorted(missing))
            )


        for line_number, row in enumerate(reader, start=2):
            
            phase_raw = clean(row["phase"])
            run_number = clean(row["run_number"])


            # ------------------------------------------------
            # Normalize phase name
            #
            # Input examples:
            #
            #   before_annealing
            #   annealing_T=75_h=5
            #   annealing_T=75_h=25
            #   annealing_T=100_h=5
            #
            # Output convention:
            #
            #   before_annealing
            #   annealing_75_5
            #   annealing_75_25
            #   annealing_100_5
            # ------------------------------------------------
            if phase_raw == "before_annealing":

                phase = "before_annealing"

            elif phase_raw.startswith("annealing_T="):

                phase = (
                    phase_raw
                    .replace("annealing_T=", "annealing_", 1)
                    .replace("_h=", "_")
                )

            else:

                raise RuntimeError(
                    f"Unrecognized phase format '{phase_raw}' "
                    f"in {filename}, line {line_number}"
                )
            if not phase:

                raise RuntimeError(
                    f"Empty phase in {filename}, "
                    f"line {line_number}"
                )


            # ------------------------------------------------
            # Construct normalized row
            # ------------------------------------------------
            normalized = {
                "sensor": sensor,
                "phase": phase,
                "run_number": run_number,
            }


            for input_column, output_column in COLUMN_MAP.items():

                normalized[output_column] = parse_integer(
                    row[input_column],
                    filename,
                    line_number,
                    input_column
                )


            key = (sensor, phase)


            if key in data:

                raise RuntimeError(
                    f"Duplicate sensor/phase in {filename}: "
                    f"{sensor}, {phase}"
                )


            data[key] = normalized
            rows.append(normalized)


    return {
        "rows": rows,
        "data": data,
    }


# ============================================================
# Read all datasets
# ============================================================

a1_24 = read_phase_file(A1_24, "A1")
a1_30 = read_phase_file(A1_30, "A1")
a1_36 = read_phase_file(A1_36, "A1")

b1_24 = read_phase_file(B1_24, "B1")
b1_30 = read_phase_file(B1_30, "B1")
b1_36 = read_phase_file(B1_36, "B1")


# ============================================================
# Print input summary
# ============================================================

print("Input rows:")
print()

print(f"A1 pixel=24 : {len(a1_24['rows'])}")
print(f"A1 pixel=30 : {len(a1_30['rows'])}")
print(f"A1 pixel=36 : {len(a1_36['rows'])}")

print()

print(f"B1 pixel=24 : {len(b1_24['rows'])}")
print(f"B1 pixel=30 : {len(b1_30['rows'])}")
print(f"B1 pixel=36 : {len(b1_36['rows'])}")

print()


# ============================================================
# Check phase correspondence
# ============================================================

def check_phases(sensor, d24, d30, d36):

    keys24 = set(d24["data"])
    keys30 = set(d30["data"])
    keys36 = set(d36["data"])

    errors = []


    for key in keys30:

        _, phase = key

        if key not in keys24:

            errors.append(
                f"{sensor}: phase '{phase}' "
                f"missing in pixel=24"
            )

        if key not in keys36:

            errors.append(
                f"{sensor}: phase '{phase}' "
                f"missing in pixel=36"
            )


    for key in keys24 - keys30:

        errors.append(
            f"{sensor}: phase '{key[1]}' exists "
            f"in pixel=24 but not in pixel=30"
        )


    for key in keys36 - keys30:

        errors.append(
            f"{sensor}: phase '{key[1]}' exists "
            f"in pixel=36 but not in pixel=30"
        )


    if errors:

        print()
        print("ERROR: phase mismatch between datasets:")
        print()

        for error in errors:
            print("  " + error)

        sys.exit(1)


check_phases("A1", a1_24, a1_30, a1_36)
check_phases("B1", b1_24, b1_30, b1_36)


# ============================================================
# Quantities for which systematic uncertainty is calculated
# ============================================================

QUANTITIES = [
    "counted_hotspots",
    "reference",
    "appeared",
    "reappeared",
    "disappeared",
    "absent",
    "not_yet_detected",
]


# ============================================================
# Systematic uncertainty
#
# delta1 = |X36 - X30|
# delta2 = |X30 - X24|
# error  = max(delta1, delta2)
# ============================================================

def systematic_error(x24, x30, x36):

    delta1 = abs(x36 - x30)
    delta2 = abs(x30 - x24)

    return max(delta1, delta2)


# ============================================================
# Process one sensor
# ============================================================

def process_sensor(sensor, d24, d30, d36):

    output_rows = []

    print(
        "============================================================"
    )
    print(f"Sensor {sensor}")
    print(
        "============================================================"
    )


    # Preserve pixel=30 phase order
    for row30 in d30["rows"]:

        phase = row30["phase"]

        key = (sensor, phase)

        row24 = d24["data"][key]
        row36 = d36["data"][key]


        # ----------------------------------------------------
        # Complete output row starts from nominal pixel=30
        # ----------------------------------------------------
        output = dict(row30)


        print()
        print(f"{sensor} - {phase}")


        # ----------------------------------------------------
        # Calculate error for all quantities
        # ----------------------------------------------------
        for quantity in QUANTITIES:

            x24 = row24[quantity]
            x30 = row30[quantity]
            x36 = row36[quantity]


            error = systematic_error(
                x24,
                x30,
                x36
            )


            output[f"{quantity}_error"] = error


            print(
                f"  {quantity:20s} "
                f"N24={x24:4d}  "
                f"N30={x30:4d}  "
                f"N36={x36:4d}  "
                f"error={error:4d}"
            )


        output_rows.append(output)


    print()

    return output_rows


# ============================================================
# Process both sensors
# ============================================================

rows_A1 = process_sensor(
    "A1",
    a1_24,
    a1_30,
    a1_36
)

rows_B1 = process_sensor(
    "B1",
    b1_24,
    b1_30,
    b1_36
)


all_rows = rows_A1 + rows_B1


# ============================================================
# OUTPUT 1
#
# Only counted hotspots
# ============================================================

summary_fieldnames = [
    "sensor",
    "phase",
    "counted_hotspots",
    "counted_hotspots_error",
]


with open(
    output_summary,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=summary_fieldnames
    )

    writer.writeheader()

    for row in all_rows:

        writer.writerow({
            field: row[field]
            for field in summary_fieldnames
        })


# ============================================================
# OUTPUT 2
#
# Complete phase-change information
# ============================================================

complete_fieldnames = [
    "sensor",
    "phase",
    "run_number",

    "counted_hotspots",
    "counted_hotspots_error",

    "reference",
    "reference_error",

    "appeared",
    "appeared_error",

    "reappeared",
    "reappeared_error",

    "disappeared",
    "disappeared_error",

    "absent",
    "absent_error",

    "not_yet_detected",
    "not_yet_detected_error",
]


with open(
    output_complete,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=complete_fieldnames
    )

    writer.writeheader()

    for row in all_rows:

        writer.writerow({
            field: row[field]
            for field in complete_fieldnames
        })


# ============================================================
# Final summary
# ============================================================

print()
print(
    "============================================================"
)
print("Files successfully written")
print(
    "============================================================"
)
print()

print("Summary:")
print(f"  {output_summary}")

print()

print("Complete:")
print(f"  {output_complete}")

print()

print(f"Total rows written: {len(all_rows)}")

PY


echo
echo "Done."