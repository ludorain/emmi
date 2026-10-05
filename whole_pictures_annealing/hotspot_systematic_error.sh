#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# 8hotspot_systematic_error.sh
#
# Calculate the systematic uncertainty on counted_hotspots
# comparing pixel sizes 24, 30 and 36.
#
# Reference:
#   N = N30
#
# Systematic uncertainty:
#   deltaN1 = |N36 - N30|
#   deltaN2 = |N30 - N24|
#   deltaN  = max(deltaN1, deltaN2)
#
# Rows are matched using:
#   sensor + phase
#
# Output:
#   hotspot_phase_summary_all_sensor_error.csv
# ============================================================


# ------------------------------------------------------------
# Directory containing this script
# ------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
FILE_24="${SCRIPT_DIR}/pixels=24/hotspot_phase_summary_all_sensors_pixel=24.csv"
FILE_30="${SCRIPT_DIR}/pixels=30/hotspot_phase_summary_all_sensors_pixel=30.csv"
FILE_36="${SCRIPT_DIR}/pixels=36/hotspot_phase_summary_all_sensors_pixel=36.csv"


# ------------------------------------------------------------
# Output file
# ------------------------------------------------------------
OUTPUT="${SCRIPT_DIR}/hotspot_phase_summary_all_sensor_error.csv"


# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------
for file in "$FILE_24" "$FILE_30" "$FILE_36"; do

    if [[ ! -f "$file" ]]; then
        echo "ERROR: input file not found:"
        echo "  $file"
        exit 1
    fi

done


echo "============================================================"
echo " Hotspot systematic uncertainty"
echo "============================================================"
echo
echo "Reference:"
echo "  pixel = 30"
echo
echo "Input files:"
echo "  $FILE_24"
echo "  $FILE_30"
echo "  $FILE_36"
echo


# ------------------------------------------------------------
# Process CSV files
# ------------------------------------------------------------
python3 - "$FILE_24" "$FILE_30" "$FILE_36" "$OUTPUT" <<'PY'

import csv
import sys
import os


file24 = sys.argv[1]
file30 = sys.argv[2]
file36 = sys.argv[3]
output = sys.argv[4]


# ============================================================
# Read one CSV and construct a dictionary indexed by
#
#     (sensor, phase)
#
# containing counted_hotspots.
# ============================================================
def read_file(filename):

    data = {}

    with open(filename, "r", newline="", encoding="utf-8-sig") as f:

        reader = csv.DictReader(f)

        required_columns = {
            "sensor",
            "phase",
            "counted_hotspots"
        }

        if reader.fieldnames is None:
            raise RuntimeError(
                f"CSV has no header: {filename}"
            )

        # Remove possible whitespace from header names
        reader.fieldnames = [
            name.strip() if name is not None else name
            for name in reader.fieldnames
        ]

        missing = required_columns - set(reader.fieldnames)

        if missing:
            raise RuntimeError(
                f"Missing columns in {filename}: "
                + ", ".join(sorted(missing))
            )


        for line_number, row in enumerate(reader, start=2):

            sensor = row["sensor"].strip()
            phase = row["phase"].strip()

            value_string = row["counted_hotspots"].strip()


            # ----------------------------------------------
            # Basic checks
            # ----------------------------------------------
            if not sensor:
                raise RuntimeError(
                    f"Empty sensor in {filename}, line {line_number}"
                )

            if not phase:
                raise RuntimeError(
                    f"Empty phase in {filename}, line {line_number}"
                )

            try:
                counted_hotspots = int(value_string)
            except ValueError:
                raise RuntimeError(
                    f"Invalid counted_hotspots='{value_string}' "
                    f"in {filename}, line {line_number}"
                )


            key = (sensor, phase)


            # ----------------------------------------------
            # Check for duplicated sensor + phase
            # ----------------------------------------------
            if key in data:
                raise RuntimeError(
                    f"Duplicate sensor/phase in {filename}: "
                    f"sensor={sensor}, phase={phase}"
                )


            data[key] = counted_hotspots


    return data


# ============================================================
# Read all three files
# ============================================================
N24 = read_file(file24)
N30 = read_file(file30)
N36 = read_file(file36)


print(f"Rows pixel=24 : {len(N24)}")
print(f"Rows pixel=30 : {len(N30)}")
print(f"Rows pixel=36 : {len(N36)}")
print()


# ============================================================
# Use pixel=30 as reference.
#
# Check that every sensor+phase present in pixel=30 also exists
# in pixel=24 and pixel=36.
# ============================================================
missing_entries = []

for key in N30:

    sensor, phase = key

    if key not in N24:
        missing_entries.append(
            f"Missing in pixel=24: sensor={sensor}, phase={phase}"
        )

    if key not in N36:
        missing_entries.append(
            f"Missing in pixel=36: sensor={sensor}, phase={phase}"
        )


if missing_entries:

    print("ERROR: sensor/phase mismatch between files:")
    print()

    for message in missing_entries:
        print("  " + message)

    sys.exit(1)


# ============================================================
# Calculate systematic uncertainty
# ============================================================
rows = []

for key, n30 in N30.items():

    sensor, phase = key

    n24 = N24[key]
    n36 = N36[key]

    deltaN1 = abs(n36 - n30)
    deltaN2 = abs(n30 - n24)

    deltaN = max(deltaN1, deltaN2)

    rows.append({
        "sensor": sensor,
        "phase": phase,
        "counted_hotspots": n30,
        "counted_hotspots_error": deltaN
    })


# ============================================================
# Write output
# ============================================================
with open(output, "w", newline="", encoding="utf-8") as f:

    fieldnames = [
        "sensor",
        "phase",
        "counted_hotspots",
        "counted_hotspots_error"
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(rows)


# ============================================================
# Print a small summary
# ============================================================
print("Systematic uncertainties:")
print()

for row in rows:

    sensor = row["sensor"]
    phase = row["phase"]
    n30 = row["counted_hotspots"]
    error = row["counted_hotspots_error"]

    n24 = N24[(sensor, phase)]
    n36 = N36[(sensor, phase)]

    print(
        f"{sensor:2s}  "
        f"{phase:25s}  "
        f"N24={n24:4d}  "
        f"N30={n30:4d}  "
        f"N36={n36:4d}  "
        f"deltaN={error:4d}"
    )


print()
print(f"Output written to:")
print(f"  {output}")

PY


echo
echo "Done."