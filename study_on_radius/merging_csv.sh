#!/bin/bash

set -euo pipefail


# ============================================================
# USER CONFIGURATION
# ============================================================

BASE_DIR="$(pwd)"

# Output CSV
OUTPUT_FILE="$BASE_DIR/DATA/bef_ann_A1_T=20_radius_study.csv"

# Maximum Euclidean distance for considering two positions
# as the same spot.
COORDINATE_TOLERANCE=5.0

# Add or remove CSV files freely.
CSV_FILES=(
    "$BASE_DIR/DATA/bef_ann_A1_T=20_r=8/luminosity/A1_T=20_total_r=8.csv"
    "$BASE_DIR/DATA/bef_ann_A1_T=20_r=10/luminosity/A1_T=20_total_r=10.csv"
    "$BASE_DIR/DATA/bef_ann_A1_T=20_r=12/luminosity/A1_T=20_total_r=12.csv"
    "$BASE_DIR/DATA/bef_ann_A1_T=20_r=14/luminosity/A1_T=20_total_r=14.csv"
    "$BASE_DIR/DATA/bef_ann_A1_T=20_r=16/luminosity/A1_T=20_total_r=16.csv"
    "$BASE_DIR/DATA/bef_ann_A1_T=20_r=18/luminosity/A1_T=20_total_r=18.csv"
    "$BASE_DIR/DATA/bef_ann_A1_T=20_r=20/luminosity/A1_T=20_total_r=20.csv"
    "$BASE_DIR/DATA/bef_ann_A1_T=20_r=22/luminosity/A1_T=20_total_r=22.csv"
    "$BASE_DIR/DATA/bef_ann_A1_T=20_r=24/luminosity/A1_T=20_total_r=24.csv"
    "$BASE_DIR/DATA/bef_ann_A1_T=20_r=26/luminosity/A1_T=20_total_r=26.csv"
)


# ============================================================
# BASH CHECKS
# ============================================================

if [[ ${#CSV_FILES[@]} -eq 0 ]]; then
    echo "Error: no input CSV files were specified."
    exit 1
fi

for csv_file in "${CSV_FILES[@]}"; do
    if [[ ! -f "$csv_file" ]]; then
        echo "Error: input CSV file not found:"
        echo "  $csv_file"
        exit 1
    fi
done

mkdir -p "$(dirname "$OUTPUT_FILE")"


# ============================================================
# EMBEDDED PYTHON CODE
# ============================================================

python3 - \
    "$OUTPUT_FILE" \
    "$COORDINATE_TOLERANCE" \
    "${CSV_FILES[@]}" <<'PYTHON_CODE'

from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    print(
        "Error: pandas is not installed.\n"
        "Install it with:\n"
        "  python3 -m pip install pandas",
        file=sys.stderr,
    )
    sys.exit(1)


# ============================================================
# INPUT ARGUMENTS FROM BASH
# ============================================================

if len(sys.argv) < 4:
    print(
        "Error: output file, tolerance and at least one "
        "input CSV are required.",
        file=sys.stderr,
    )
    sys.exit(1)

output_file = Path(sys.argv[1]).expanduser().resolve()
coordinate_tolerance = float(sys.argv[2])

input_files = [
    Path(filename).expanduser().resolve()
    for filename in sys.argv[3:]
]


REQUIRED_COLUMNS = [
    "spot",
    "x",
    "y",
    "luminosity",
    "error",
    "T",
    "v",
    "v_fin",
    "integration_radius",
]


# ============================================================
# UNION-FIND
# ============================================================

class UnionFind:
    """
    Data structure used to join positions belonging to the
    same global spot.
    """

    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, index: int) -> int:
        if self.parent[index] != index:
            self.parent[index] = self.find(self.parent[index])

        return self.parent[index]

    def union(self, first: int, second: int) -> None:
        root_first = self.find(first)
        root_second = self.find(second)

        if root_first == root_second:
            return

        if self.rank[root_first] < self.rank[root_second]:
            self.parent[root_first] = root_second

        elif self.rank[root_first] > self.rank[root_second]:
            self.parent[root_second] = root_first

        else:
            self.parent[root_second] = root_first
            self.rank[root_first] += 1


# ============================================================
# READ THE CSV FILES
# ============================================================

if coordinate_tolerance <= 0:
    print(
        "Error: COORDINATE_TOLERANCE must be greater than zero.",
        file=sys.stderr,
    )
    sys.exit(1)

dataframes = []

for source_index, input_file in enumerate(input_files):

    print(f"Reading: {input_file}")

    try:
        dataframe = pd.read_csv(input_file)

    except Exception as error:
        print(
            f"Error while reading '{input_file}': {error}",
            file=sys.stderr,
        )
        sys.exit(1)

    missing_columns = [
        column
        for column in REQUIRED_COLUMNS
        if column not in dataframe.columns
    ]

    if missing_columns:
        print(
            f"Error: file '{input_file}' is missing columns:\n"
            f"  {', '.join(missing_columns)}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Convert the required columns to numeric values.
    for column in REQUIRED_COLUMNS:
        dataframe[column] = pd.to_numeric(
            dataframe[column],
            errors="coerce",
        )

    invalid_rows = dataframe[REQUIRED_COLUMNS].isna().any(axis=1)

    if invalid_rows.any():
        csv_row_numbers = (
            dataframe.index[invalid_rows] + 2
        ).tolist()

        print(
            f"Error: non-numeric or missing values in "
            f"'{input_file}'.\n"
            f"Problematic CSV rows: {csv_row_numbers}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Save the original local spot and file index.
    dataframe["_source_index"] = source_index
    dataframe["_original_spot"] = dataframe["spot"].astype(int)

    dataframes.append(dataframe)


merged = pd.concat(
    dataframes,
    ignore_index=True,
    sort=False,
)


# ============================================================
# REPRESENTATIVE POSITION OF EACH LOCAL SPOT
# ============================================================

# Each original spot normally appears several times, once for
# every v_fin or integration_radius value.
#
# Therefore, first create one representative x-y position for
# each spot in each input file.

local_spots = (
    merged
    .groupby(
        ["_source_index", "_original_spot"],
        as_index=False,
        sort=False,
    )
    .agg(
        representative_x=("x", "median"),
        representative_y=("y", "median"),
    )
)


number_of_local_spots = len(local_spots)

if number_of_local_spots == 0:
    print(
        "Error: no spots were found in the input files.",
        file=sys.stderr,
    )
    sys.exit(1)


# ============================================================
# MATCH SPOTS USING THE X-Y DISTANCE
# ============================================================

union_find = UnionFind(number_of_local_spots)

x_values = local_spots["representative_x"].to_numpy(dtype=float)
y_values = local_spots["representative_y"].to_numpy(dtype=float)

# Spatial grid used to reduce the number of distance comparisons.
spatial_grid: dict[tuple[int, int], list[int]] = defaultdict(list)

for current_index in range(number_of_local_spots):

    x_current = x_values[current_index]
    y_current = y_values[current_index]

    grid_x = math.floor(x_current / coordinate_tolerance)
    grid_y = math.floor(y_current / coordinate_tolerance)

    # Only points in the same grid cell or in adjacent cells
    # can be separated by less than the tolerance.
    for delta_x in (-1, 0, 1):
        for delta_y in (-1, 0, 1):

            neighbour_cell = (
                grid_x + delta_x,
                grid_y + delta_y,
            )

            for other_index in spatial_grid.get(
                neighbour_cell,
                [],
            ):

                distance = math.hypot(
                    x_current - x_values[other_index],
                    y_current - y_values[other_index],
                )

                if distance <= coordinate_tolerance:
                    union_find.union(
                        current_index,
                        other_index,
                    )

    spatial_grid[(grid_x, grid_y)].append(current_index)


# ============================================================
# CREATE THE NEW GLOBAL SPOT INDICES
# ============================================================

components: dict[int, list[int]] = defaultdict(list)

for local_index in range(number_of_local_spots):
    component_root = union_find.find(local_index)
    components[component_root].append(local_index)


# Calculate the average position of every global spot.
component_centers = []

for component_root, indices in components.items():

    mean_x = sum(
        x_values[index]
        for index in indices
    ) / len(indices)

    mean_y = sum(
        y_values[index]
        for index in indices
    ) / len(indices)

    component_centers.append(
        (
            component_root,
            mean_x,
            mean_y,
        )
    )


# Deterministic numbering:
# first by increasing x coordinate, then by increasing y.
component_centers.sort(
    key=lambda element: (
        element[1],
        element[2],
    )
)

component_to_global_spot = {
    component_root: global_spot
    for global_spot, (
        component_root,
        _,
        _,
    ) in enumerate(component_centers)
}


# Connect each local spot to its new global spot index.
local_to_global: dict[tuple[int, int], int] = {}

for local_index, row in local_spots.iterrows():

    component_root = union_find.find(local_index)

    global_spot = component_to_global_spot[
        component_root
    ]

    key = (
        int(row["_source_index"]),
        int(row["_original_spot"]),
    )

    local_to_global[key] = global_spot


# Assign the new spot value to all rows.
new_spot_values = []

for source_index, original_spot in zip(
    merged["_source_index"],
    merged["_original_spot"],
):

    key = (
        int(source_index),
        int(original_spot),
    )

    new_spot_values.append(
        local_to_global[key]
    )

merged["spot"] = new_spot_values
merged["spot"] = merged["spot"].astype(int)


# ============================================================
# SORT THE OUTPUT
# ============================================================

merged = merged.sort_values(
    by=[
        "spot",
        "v_fin",
        "integration_radius",
    ],
    ascending=[
        True,
        True,
        True,
    ],
    kind="mergesort",
)


# Remove internal helper columns.
merged = merged.drop(
    columns=[
        "_source_index",
        "_original_spot",
    ]
)


# Put the standard columns first.
first_columns = [
    "spot",
    "x",
    "y",
    "luminosity",
    "error",
    "T",
    "v",
    "v_fin",
    "integration_radius",
]

additional_columns = [
    column
    for column in merged.columns
    if column not in first_columns
]

merged = merged[
    first_columns + additional_columns
]

merged = merged.reset_index(drop=True)


# ============================================================
# SAVE THE OUTPUT CSV
# ============================================================

output_file.parent.mkdir(
    parents=True,
    exist_ok=True,
)

merged.to_csv(
    output_file,
    index=False,
)

print()
print("Merge completed successfully.")
print(f"Number of input CSV files: {len(input_files)}")
print(f"Total number of rows: {len(merged)}")
print(f"Number of global spots: {merged['spot'].nunique()}")
print(f"Coordinate tolerance: {coordinate_tolerance}")
print(f"Output file: {output_file}")

PYTHON_CODE


echo
echo "The merged CSV was saved in:"
echo "$OUTPUT_FILE"