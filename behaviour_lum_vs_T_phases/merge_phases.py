#!/usr/bin/env python3

from pathlib import Path

import pandas as pd


# ============================================================
# INPUT AND OUTPUT PATHS
# ============================================================

BEFORE_ANNEALING_CSV = Path(
    "/Users/ludovicarainero/emmi/DATA_irradiated/merged_files_v=const/A1_v=5_run=20260514-062228_total_global_ID.csv"
)

ANNEALING_75_5_CSV = Path(
    "/path/to/annealing_T=75_h=5.csv"
)

ANNEALING_75_25_CSV = Path(
    "/path/to/annealing_T=75_h=25.csv"
)

ANNEALING_100_5_CSV = Path(
    "/path/to/annealing_T=100_h=5.csv"
)

ANNEALING_100_25_CSV = Path(
    "/path/to/annealing_T=100_h=25.csv"
)

OUTPUT_CSV = Path(
    "/Users/ludovicarainero/emmi/behaviour_lum_vs_T_phases/lum_vs_T_complete.csv"
)


# ============================================================
# CONFIGURATION
# ============================================================

REQUIRED_COLUMNS = [
    "spot",
    "x",
    "y",
    "luminosity",
    "error",
    "T",
    "v",
]

# The order used when sorting the phase column
PHASES = [
    "before_annealing",
    "annealing_T=75_h=5",
    "annealing_T=75_h=25",
    "annealing_T=100_h=5",
    "annealing_T=100_h=25",
]

INPUT_FILES = [
    (BEFORE_ANNEALING_CSV, "before_annealing"),
    (ANNEALING_75_5_CSV, "annealing_T=75_h=5"),
    (ANNEALING_75_25_CSV, "annealing_T=75_h=25"),
    (ANNEALING_100_5_CSV, "annealing_T=100_h=5"),
    (ANNEALING_100_25_CSV, "annealing_T=100_h=25"),
]


def read_and_add_phase(csv_path, phase):
    """Read one CSV file, check its columns, and add the phase."""

    if not csv_path.is_file():
        raise FileNotFoundError(
            f"Input file not found:\n{csv_path}"
        )

    dataframe = pd.read_csv(csv_path)

    missing_columns = [
        column
        for column in REQUIRED_COLUMNS
        if column not in dataframe.columns
    ]

    if missing_columns:
        raise ValueError(
            f"File '{csv_path}' is missing the columns: "
            f"{', '.join(missing_columns)}"
        )

    # Keep only the expected columns
    dataframe = dataframe[REQUIRED_COLUMNS].copy()

    # Add the phase column
    dataframe["phase"] = phase

    return dataframe


def main():

    dataframes = []

    for csv_path, phase in INPUT_FILES:
        dataframe = read_and_add_phase(
            csv_path=csv_path,
            phase=phase,
        )

        dataframes.append(dataframe)

        print(
            f"Loaded {csv_path.name}: "
            f"{len(dataframe)} rows, phase = {phase}"
        )

    # Merge the five CSV files
    merged_dataframe = pd.concat(
        dataframes,
        ignore_index=True,
    )

    # Define the chronological order of the phases
    merged_dataframe["phase"] = pd.Categorical(
        merged_dataframe["phase"],
        categories=PHASES,
        ordered=True,
    )

    # Sort by spot, temperature, and phase
    merged_dataframe = merged_dataframe.sort_values(
        by=["spot", "T", "phase"],
        ascending=[True, True, True],
        kind="stable",
    ).reset_index(drop=True)

    # Convert phase back to normal strings
    merged_dataframe["phase"] = merged_dataframe[
        "phase"
    ].astype(str)

    # Create the output directory if it does not exist
    OUTPUT_CSV.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Save the merged CSV
    merged_dataframe.to_csv(
        OUTPUT_CSV,
        index=False,
    )

    print()
    print(f"Merged CSV saved to:\n{OUTPUT_CSV}")
    print(f"Total number of rows: {len(merged_dataframe)}")


if __name__ == "__main__":
    main()