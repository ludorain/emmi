
#!/usr/bin/env bash

# ============================================================
# STUDY OF LUMINOSITY VS INTEGRATION RADIUS
#
# Usage:
#
#   ./study_luminosity_vs_radius.sh A1 T before_annealing
#   ./study_luminosity_vs_radius.sh A1 v annealing_T=75_h=5
#
# The script must be located in:
#
#   emmi/study_on_radius/
#
# For every selected hotspot it:
#
#   1. reads its GLOBAL ID and x,y coordinates from the
#      geometry plan produced by automatic_irradiated_analysis.sh
#
#   2. ignores the original integration radius
#
#   3. scans:
#
#          r = 2, 3, 4, ..., 25 pixels
#
#   4. creates a ROOT coordinate file with
#
#          x, y, pi*r^2, r
#
#   5. runs spot_luminosity_sum_irradiated.C
#      on every TH2F file of the selected phase
#
#   6. creates one CSV for every radius
#
#   7. merges all radii into one final CSV
#
# Final columns:
#
#   ID_hotspot_globale
#   x
#   y
#   integration_radius
#   luminosity
#   error
#   T
#   v
#
# ============================================================

set -Eeuo pipefail
shopt -s nullglob


# ============================================================
# BASIC FUNCTIONS
# ============================================================

die() {
    echo "ERROR: $*" >&2
    exit 1
}


log() {
    echo
    echo "============================================================"
    echo "$*"
    echo "============================================================"
}


# ============================================================
# INPUT
# ============================================================

if [[ $# -ne 3 ]]; then

    cat >&2 <<EOF

Usage:

    $0 SENSOR CONSTANT PHASE

Examples:

    $0 A1 T before_annealing
    $0 A1 T annealing_T=75_h=5
    $0 A1 v before_annealing

SENSOR:
    A1 | A2 | B1 | B2

CONSTANT:
    T | v

PHASE:
    before_annealing
    annealing_T=75_h=5
    annealing_T=75_h=25
    annealing_T=100_h=5
    annealing_T=100_h=25
    ...

EOF

    exit 1
fi


SENSOR="$1"
CONSTANT="$2"
PHASE="$3"


case "$SENSOR" in
    A1|A2|B1|B2)
        ;;
    *)
        die "Invalid sensor '$SENSOR'."
        ;;
esac


case "$CONSTANT" in
    T|v)
        ;;
    *)
        die "Invalid constant '$CONSTANT'. Use T or v."
        ;;
esac


# ============================================================
# PATHS
# ============================================================

# This script is inside emmi/study_on_radius
STUDY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Main emmi directory
EMMI_DIR="$(cd "$STUDY_DIR/.." && pwd)"


# ------------------------------------------------------------
# Results produced by yesterday's automatic analysis
# ------------------------------------------------------------

MERGED_ROOT="$EMMI_DIR/DATA_irradiated_isolated_changeR/merged_files"


# ------------------------------------------------------------
# TH2F data that must be used for the radius study
# ------------------------------------------------------------

TH2F_DATA_ROOT="$EMMI_DIR/DATA_irradiated_isolated_changeR"


# ------------------------------------------------------------
# ROOT luminosity code
# ------------------------------------------------------------

SPOT_LUM_DIR="$EMMI_DIR/spot_luminosity"

SPOT_LUM_MACRO="spot_luminosity_sum_irradiated.C"


# ------------------------------------------------------------
# Python
# ------------------------------------------------------------

PYTHON="${PYTHON:-python3}"


# ============================================================
# CHECK REQUIREMENTS
# ============================================================

command -v "$PYTHON" >/dev/null 2>&1 \
    || die "Python not found."

command -v root >/dev/null 2>&1 \
    || die "ROOT not found."

[[ -d "$MERGED_ROOT" ]] \
    || die "Merged-files directory not found: $MERGED_ROOT"

[[ -d "$TH2F_DATA_ROOT" ]] \
    || die "TH2F data directory not found: $TH2F_DATA_ROOT"

[[ -f "$SPOT_LUM_DIR/$SPOT_LUM_MACRO" ]] \
    || die "ROOT macro not found: $SPOT_LUM_DIR/$SPOT_LUM_MACRO"


# ============================================================
# FIND SENSOR/CONDITION DIRECTORY
#
# Example:
#
# DATA_irradiated/merged_files/A1_T=20
#
# We deliberately do NOT hard-code T=20 or v=5.
# The directory produced by yesterday's code determines it.
# ============================================================

MERGED_MATCHES=(
    "$MERGED_ROOT/${SENSOR}_${CONSTANT}="*
)


if (( ${#MERGED_MATCHES[@]} == 0 )); then

    die \
        "No merged directory found for ${SENSOR}_${CONSTANT}=* in:
$MERGED_ROOT"

fi


if (( ${#MERGED_MATCHES[@]} > 1 )); then

    echo "More than one merged directory found:" >&2

    printf '  %s\n' "${MERGED_MATCHES[@]}" >&2

    die "Sensor + constant must identify exactly one dataset."

fi


MERGED_DIR="${MERGED_MATCHES[0]}"

DATASET_LABEL="$(basename "$MERGED_DIR")"


# Example:
#
# DATASET_LABEL = A1_T=20


# ============================================================
# GEOMETRY PLAN
# ============================================================

GEOMETRY_PLAN="$MERGED_DIR/${SENSOR}_${CONSTANT}_geometry_plan.csv"


[[ -f "$GEOMETRY_PLAN" ]] \
    || die "Geometry plan not found: $GEOMETRY_PLAN"


# ============================================================
# FIND TH2F RUN
#
# Example:
#
# DATA_irradiated_isolated_changeR/
# └── before_annealing/
#     └── A1_T=20_run=20260513-032137/
#         └── 5th2f/
# ============================================================

PHASE_DIR="$TH2F_DATA_ROOT/$PHASE"


[[ -d "$PHASE_DIR" ]] \
    || die "Phase directory not found: $PHASE_DIR"


RUN_MATCHES=(
    "$PHASE_DIR/${DATASET_LABEL}_run="*
)


VALID_RUNS=()

for run_dir in "${RUN_MATCHES[@]}"; do

    if [[ -d "$run_dir/5th2f" ]]; then
        VALID_RUNS+=("$run_dir")
    fi

done


if (( ${#VALID_RUNS[@]} == 0 )); then

    die \
        "No run with 5th2f found for:

${DATASET_LABEL}

inside:

${PHASE_DIR}"

fi


if (( ${#VALID_RUNS[@]} > 1 )); then

    echo "More than one compatible run found:" >&2

    printf '  %s\n' "${VALID_RUNS[@]}" >&2

    die "Expected exactly one run for the selected phase."

fi


RUN_DIR="${VALID_RUNS[0]}"
TH2F_DIR="$RUN_DIR/5th2f"


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR="$STUDY_DIR/${DATASET_LABEL}_${PHASE}_luminosities"

mkdir -p "$OUTPUT_DIR"


# Temporary working directory
TMP_DIR="$OUTPUT_DIR/.radius_scan_tmp"

rm -rf "$TMP_DIR"

mkdir -p "$TMP_DIR"


# ============================================================
# SUMMARY
# ============================================================

log "LUMINOSITY VS INTEGRATION RADIUS"

echo "Sensor:          $SENSOR"
echo "Constant:        $CONSTANT"
echo "Dataset:         $DATASET_LABEL"
echo "Phase:           $PHASE"
echo
echo "Geometry plan:"
echo "  $GEOMETRY_PLAN"
echo
echo "TH2F directory:"
echo "  $TH2F_DIR"
echo
echo "Output directory:"
echo "  $OUTPUT_DIR"
echo
echo "Radius scan:"
echo "  2 -> 25 pixels"


# ============================================================
# MAIN ANALYSIS
# ============================================================

SENSOR_ENV="$SENSOR" \
CONSTANT_ENV="$CONSTANT" \
PHASE_ENV="$PHASE" \
DATASET_LABEL_ENV="$DATASET_LABEL" \
GEOMETRY_PLAN_ENV="$GEOMETRY_PLAN" \
TH2F_DIR_ENV="$TH2F_DIR" \
OUTPUT_DIR_ENV="$OUTPUT_DIR" \
TMP_DIR_ENV="$TMP_DIR" \
SPOT_LUM_DIR_ENV="$SPOT_LUM_DIR" \
SPOT_LUM_MACRO_ENV="$SPOT_LUM_MACRO" \
"$PYTHON" <<'PYCODE'

from pathlib import Path

import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd


# ============================================================
# CONFIGURATION FROM BASH
# ============================================================

sensor = os.environ["SENSOR_ENV"]

constant = os.environ["CONSTANT_ENV"]

phase = os.environ["PHASE_ENV"]

dataset_label = os.environ["DATASET_LABEL_ENV"]

geometry_plan_path = Path(
    os.environ["GEOMETRY_PLAN_ENV"]
)

th2f_dir = Path(
    os.environ["TH2F_DIR_ENV"]
)

output_dir = Path(
    os.environ["OUTPUT_DIR_ENV"]
)

tmp_dir = Path(
    os.environ["TMP_DIR_ENV"]
)

spot_lum_dir = Path(
    os.environ["SPOT_LUM_DIR_ENV"]
)

spot_lum_macro = os.environ[
    "SPOT_LUM_MACRO_ENV"
]


# ============================================================
# RADIUS RANGE
# ============================================================

RADIUS_MIN = 2

RADIUS_MAX = 25


# ============================================================
# COORDINATE MATCHING TOLERANCE
#
# ROOT returns the x,y coordinates present in the coordinate
# file. A tolerance is nevertheless used to make the mapping
# robust against small numerical-formatting differences.
# ============================================================

COORD_MATCH_RADIUS = 1.0


# ============================================================
# IMPORTANT:
#
# The geometry plan contains also propagated coordinates for
# hotspots that were NOT detected in a given phase.
#
# For the radius study we want the hotspots genuinely detected
# in the selected phase, equivalent to using that phase's
# original coordinate catalogue.
# ============================================================

ONLY_DETECTED = True


# ============================================================
# HELPERS
# ============================================================

def boolean_mask(series):

    """
    Robust conversion of detection column to boolean.
    """

    if pd.api.types.is_bool_dtype(series):
        return series

    text = (
        series
        .astype(str)
        .str.strip()
        .str.lower()
    )

    return text.isin(
        [
            "present",
            "detected",
        ]
    )


def extract_parameter(filename, parameter):

    """
    Extract T or v from a TH2F filename.

    Example:

        ..._T=20_v=58.30_data=diff_...root
    """

    pattern = re.compile(
        rf"(?:^|_){re.escape(parameter)}="
        rf"([-+]?(?:\d+(?:\.\d*)?|\.\d+)"
        rf"(?:[eE][-+]?\d+)?)"
    )

    match = pattern.search(filename)

    if match is None:

        raise RuntimeError(
            f"Cannot extract {parameter}=... "
            f"from filename:\n{filename}"
        )

    return float(match.group(1))


def match_root_results(
    result,
    coordinates,
    root_file,
):

    """
    Associate every ROOT result with its GLOBAL hotspot.

    Matching is spatial and one-to-one.

    This prevents a dependence on the local ordering used
    internally by ROOT.
    """

    required = [
        "x",
        "y",
        "luminosity",
        "error",
    ]

    missing = [
        c
        for c in required
        if c not in result.columns
    ]

    if missing:

        raise RuntimeError(
            f"ROOT output for {root_file.name} "
            f"is missing columns {missing}.\n"
            f"Found columns: {list(result.columns)}"
        )


    result = result.copy()

    for column in required:
        result[column] = pd.to_numeric(
            result[column],
            errors="raise",
        )


    # --------------------------------------------------------
    # Number of rows must match number of requested hotspots
    # --------------------------------------------------------

    if len(result) != len(coordinates):

        raise RuntimeError(
            "\nNumber of ROOT results does not match "
            "the number of coordinates.\n\n"
            f"ROOT file: {root_file}\n"
            f"Coordinates: {len(coordinates)}\n"
            f"ROOT results: {len(result)}\n"
        )


    # --------------------------------------------------------
    # Build every possible pair inside tolerance
    # --------------------------------------------------------

    pairs = []

    for result_index, result_row in result.iterrows():

        for coord_index, coord_row in coordinates.iterrows():

            distance = math.hypot(
                result_row["x"] - coord_row["x"],
                result_row["y"] - coord_row["y"],
            )

            if distance <= COORD_MATCH_RADIUS:

                pairs.append(
                    (
                        distance,
                        result_index,
                        coord_index,
                    )
                )


    # Closest pairs first
    pairs.sort(
        key=lambda item: item[0]
    )


    assigned_results = set()

    assigned_coordinates = set()

    mapping = {}


    for distance, result_index, coord_index in pairs:

        if result_index in assigned_results:
            continue

        if coord_index in assigned_coordinates:
            continue

        mapping[result_index] = coord_index

        assigned_results.add(
            result_index
        )

        assigned_coordinates.add(
            coord_index
        )


    # --------------------------------------------------------
    # Every ROOT row must have exactly one global hotspot
    # --------------------------------------------------------

    if len(mapping) != len(coordinates):

        raise RuntimeError(
            "\nCould not associate every ROOT result "
            "with a global hotspot.\n\n"
            f"ROOT file: {root_file}\n"
            f"Matched: {len(mapping)} / {len(coordinates)}\n"
            f"Tolerance: {COORD_MATCH_RADIUS} pixels\n"
        )


    # --------------------------------------------------------
    # Construct final table
    # --------------------------------------------------------

    records = []

    for result_index, coord_index in mapping.items():

        r = result.loc[result_index]

        c = coordinates.loc[coord_index]

        records.append(
            {
                "ID_hotspot_globale":
                    int(c["ID_hotspot_globale"]),

                # Use authoritative coordinates from
                # the geometry plan.
                "x":
                    float(c["x"]),

                "y":
                    float(c["y"]),

                "luminosity":
                    float(r["luminosity"]),

                "error":
                    float(r["error"]),
            }
        )


    return pd.DataFrame(records)


# ============================================================
# READ GEOMETRY PLAN
# ============================================================

print()
print("Reading geometry plan:")
print(geometry_plan_path)


geometry = pd.read_csv(
    geometry_plan_path
)


required_columns = [
    "spot",
    "phase",
    "x",
    "y",
    "integration_radius",
]


missing = [
    c
    for c in required_columns
    if c not in geometry.columns
]


if missing:

    raise RuntimeError(
        "Geometry plan is missing columns: "
        + ", ".join(missing)
    )


# ============================================================
# SELECT PHASE
# ============================================================

phase_geometry = geometry[
    geometry["phase"].astype(str) == phase
].copy()


if phase_geometry.empty:

    available_phases = sorted(
        geometry["phase"]
        .astype(str)
        .unique()
    )

    raise RuntimeError(
        f"\nPhase '{phase}' was not found "
        f"in the geometry plan.\n\n"
        "Available phases:\n  "
        + "\n  ".join(available_phases)
    )


# ============================================================
# KEEP ONLY REAL DETECTIONS
# ============================================================

if ONLY_DETECTED and "detection" in phase_geometry.columns:

    detection_mask = boolean_mask(
        phase_geometry["detection"]
    )

    phase_geometry = phase_geometry[
        detection_mask
    ].copy()


if phase_geometry.empty:

    raise RuntimeError(
        f"No detected hotspots remain "
        f"for phase {phase}."
    )


# ============================================================
# CREATE GLOBAL COORDINATE TABLE
# ============================================================

coordinates = (
    phase_geometry[
        [
            "spot",
            "x",
            "y",
            "integration_radius",
        ]
    ]
    .drop_duplicates(
        subset=["spot"]
    )
    .rename(
        columns={
            "spot":
                "ID_hotspot_globale"
        }
    )
)


coordinates[
    "ID_hotspot_globale"
] = pd.to_numeric(
    coordinates[
        "ID_hotspot_globale"
    ],
    errors="raise",
).astype(int)


coordinates["x"] = pd.to_numeric(
    coordinates["x"],
    errors="raise",
)

coordinates["y"] = pd.to_numeric(
    coordinates["y"],
    errors="raise",
)

coordinates["integration_radius"] = pd.to_numeric(
    coordinates["integration_radius"],
    errors="raise",
)


coordinates = (
    coordinates
    .sort_values(
        "ID_hotspot_globale"
    )
    .reset_index(drop=True)
)


print()
print(
    f"Selected phase: {phase}"
)

print(
    f"Hotspots used for radius study: "
    f"{len(coordinates)}"
)


# ============================================================
# SAVE REFERENCE GLOBAL COORDINATES
#
# Useful for checking exactly which spots were analysed.
# ============================================================

reference_coordinates_csv = (
    output_dir
    / f"{dataset_label}_{phase}_"
      "reference_global_coordinates.csv"
)


coordinates.to_csv(
    reference_coordinates_csv,
    index=False,
)


print(
    "Reference coordinate table:"
)

print(
    reference_coordinates_csv
)


# ============================================================
# FIND ROOT FILES
# ============================================================

root_files = sorted(
    th2f_dir.glob(
        "*data=diff*_processed_rotated_th2f.root"
    )
)


# Explicitly remove uncertainty images if their name happens
# to contain the data=diff substring.
root_files = [
    path
    for path in root_files
    if "data=diffe" not in path.name
]


if not root_files:

    raise RuntimeError(
        f"No TH2F ROOT files found in:\n"
        f"{th2f_dir}"
    )


print()
print(
    f"TH2F files found: {len(root_files)}"
)


for path in root_files:
    print(
        f"  {path.name}"
    )


# ============================================================
# OUTPUT ACCUMULATOR
# ============================================================

all_radius_frames = []


# ============================================================
# LOOP OVER INTEGRATION RADIUS
# ============================================================

for radius in range(
    RADIUS_MIN,
    RADIUS_MAX + 1,
):

    print()
    print(
        "============================================================"
    )

    print(
        f"INTEGRATION RADIUS = {radius} pixels"
    )

    print(
        "============================================================"
    )


    # --------------------------------------------------------
    # Create coordinate file for ROOT
    #
    # Format:
    #
    # x, y, integration_area, integration_radius
    # --------------------------------------------------------

    integration_area = (
        math.pi
        * radius
        * radius
    )


    coordinate_file = (
        tmp_dir
        / (
            f"{dataset_label}_"
            f"{phase}_"
            f"r={radius}_coordinates.txt"
        )
    )


    with coordinate_file.open(
        "w",
        encoding="utf-8",
    ) as handle:

        for row in coordinates.itertuples(
            index=False
        ):

            handle.write(
                f"{row.x:.6f}, "
                f"{row.y:.6f}, "
                f"{integration_area:.6f}, "
                f"{float(radius):.6f}\n"
            )


    # --------------------------------------------------------
    # Measurements from every operating condition
    # --------------------------------------------------------

    radius_measurements = []


    for root_file in root_files:

        # ----------------------------------------------------
        # Read T and v from filename
        # ----------------------------------------------------

        t_value = extract_parameter(
            root_file.name,
            "T",
        )

        v_value = extract_parameter(
            root_file.name,
            "v",
        )


        print(
            f"  TH2F: {root_file.name}"
        )

        print(
            f"        T = {t_value:g}, "
            f"v = {v_value:g}"
        )


        # ----------------------------------------------------
        # ROOT creates luminosity_results.csv inside
        # emmi/spot_luminosity
        # ----------------------------------------------------

        result_csv = (
            spot_lum_dir
            / "luminosity_results.csv"
        )


        if result_csv.exists():
            result_csv.unlink()


        macro_call = (
            f'{spot_lum_macro}'
            f'("{root_file}",'
            f'"{coordinate_file}")'
        )


        completed = subprocess.run(
            [
                "root",
                "-l",
                "-q",
                macro_call,
            ],
            cwd=spot_lum_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )


        if completed.returncode != 0:

            print(
                completed.stdout,
                file=sys.stderr,
            )

            raise RuntimeError(
                f"ROOT failed for:\n"
                f"{root_file}\n"
                f"radius = {radius}"
            )


        if not result_csv.is_file():

            print(
                completed.stdout,
                file=sys.stderr,
            )

            raise RuntimeError(
                "ROOT did not create "
                "luminosity_results.csv for:\n"
                f"{root_file}"
            )


        # ----------------------------------------------------
        # Read ROOT output
        # ----------------------------------------------------

        raw_result = pd.read_csv(
            result_csv
        )


        # Remove immediately to avoid mixing subsequent runs.
        result_csv.unlink()


        # ----------------------------------------------------
        # Recover GLOBAL hotspot IDs
        # ----------------------------------------------------

        measured = match_root_results(
            result=raw_result,
            coordinates=coordinates,
            root_file=root_file,
        )


        # ----------------------------------------------------
        # Add radius and operating conditions
        # ----------------------------------------------------

        measured[
            "integration_radius"
        ] = radius

        measured["T"] = t_value

        measured["v"] = v_value


        radius_measurements.append(
            measured
        )


    # ========================================================
    # MERGE ALL OPERATING CONDITIONS FOR THIS RADIUS
    # ========================================================

    radius_df = pd.concat(
        radius_measurements,
        ignore_index=True,
    )


    # Exact requested column order
    radius_df = radius_df[
        [
            "ID_hotspot_globale",
            "x",
            "y",
            "integration_radius",
            "luminosity",
            "error",
            "T",
            "v",
        ]
    ]


    # Good ordering for later Lum(radius) plots:
    #
    # hotspot -> operating condition -> radius
    radius_df = (
        radius_df
        .sort_values(
            [
                "ID_hotspot_globale",
                "T",
                "v",
            ],
            kind="stable",
        )
        .reset_index(drop=True)
    )


    # --------------------------------------------------------
    # Save one CSV per radius
    # --------------------------------------------------------

    radius_output = (
        output_dir
        / (
            f"{dataset_label}_"
            f"{phase}_"
            f"r={radius}.csv"
        )
    )


    radius_df.to_csv(
        radius_output,
        index=False,
    )


    print(
        f"  Created: {radius_output.name}"
    )


    all_radius_frames.append(
        radius_df
    )


# ============================================================
# FINAL MERGE: ALL RADII
# ============================================================

print()
print(
    "============================================================"
)

print(
    "CREATING FINAL ALL-RADII CSV"
)

print(
    "============================================================"
)


final_df = pd.concat(
    all_radius_frames,
    ignore_index=True,
)


# Arrange data so that, for every hotspot and operating
# condition, r goes 2 -> 25 consecutively.
final_df = (
    final_df
    .sort_values(
        [
            "ID_hotspot_globale",
            "T",
            "v",
            "integration_radius",
        ],
        kind="stable",
    )
    .reset_index(drop=True)
)


final_output = (
    output_dir
    / (
        f"{dataset_label}_"
        f"{phase}_"
        "all_radii.csv"
    )
)


final_df.to_csv(
    final_output,
    index=False,
)


# ============================================================
# CLEAN TEMPORARY FILES
# ============================================================

shutil.rmtree(
    tmp_dir,
    ignore_errors=True,
)


# ============================================================
# REPORT
# ============================================================

print()
print(
    "============================================"
)

print(
    "RADIUS STUDY COMPLETED"
)

print(
    "============================================"
)

print(
    f"Sensor:             {sensor}"
)

print(
    f"Condition:          {constant}"
)

print(
    f"Phase:              {phase}"
)

print(
    f"Number of hotspots: {len(coordinates)}"
)

print(
    f"Radius range:       "
    f"{RADIUS_MIN} - {RADIUS_MAX}"
)

print(
    f"TH2F files:         {len(root_files)}"
)

print()

print(
    "Final CSV:"
)

print(
    final_output
)

print(
    "============================================"
)

PYCODE


log "DONE"

echo "Results saved in:"
echo
echo "  $OUTPUT_DIR"
echo