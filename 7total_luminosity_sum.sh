#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 7total_luminosity_sum.sh
#
# Compute, for each available annealing phase, the mean hotspot
# luminosity at a fixed final overvoltage:
#
#   mean_luminosity = sum_i(L_i) / N_spots
#
# with statistical uncertainty:
#
#   error = sqrt(sum_i(sigma_i^2)) / N_spots
#
# Supported sensors:
#   A1 -> v_fin = 5 V
#   B1 -> v_fin = 3 V
#
# Usage:
#   ./7total_luminosity_sum.sh A1
#   ./7total_luminosity_sum.sh B1
#   ./7total_luminosity_sum.sh A1 20
#
# The temperature defaults to T = 20 C.
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SENSOR="${1:-}"
TEMPERATURE="${2:-20}"

if [[ -z "$SENSOR" ]]; then
    echo "Usage: $0 SENSOR [TEMPERATURE]"
    echo "Example: $0 A1"
    exit 1
fi

case "$SENSOR" in
    A1)
        TARGET_VFIN="5"
        ;;
    B1)
        TARGET_VFIN="3"
        ;;
    *)
        echo "ERROR: unsupported sensor '$SENSOR'."
        echo "Supported sensors: A1, B1"
        exit 1
        ;;
esac

BASE_DIR="${SCRIPT_DIR}/DATA_irradiated_isolated_R=20"
OUTPUT_DIR="${SCRIPT_DIR}/analysis/total_luminosity_sum"

if [[ ! -d "$BASE_DIR" ]]; then
    echo "ERROR: input directory not found:"
    echo "  $BASE_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

python3 - "$BASE_DIR" "$OUTPUT_DIR" "$SENSOR" "$TEMPERATURE" "$TARGET_VFIN" <<'PY'
import csv
import glob
import math
import os
import re
import sys

base_dir, output_dir, sensor, temperature_arg, target_vfin_arg = sys.argv[1:6]

try:
    temperature = float(temperature_arg)
    target_vfin = float(target_vfin_arg)
except ValueError:
    raise SystemExit("ERROR: temperature and v_fin must be numeric.")

def number_token(x):
    """Use 20 instead of 20.0 in path/output names when possible."""
    if float(x).is_integer():
        return str(int(x))
    return f"{x:g}"

temperature_token = number_token(temperature)
vfin_token = number_token(target_vfin)

def phase_sort_key(phase):
    """
    Sort phases as:
      before_annealing
      75 C, 5 h
      75 C, 25 h
      100 C, 5 h
      ...
    Unknown phase names are placed afterwards.
    """
    if phase in ("before_annealing", "bef_ann"):
        return (0, 0, 0, phase)

    m = re.fullmatch(r"annealing_T=([0-9.]+)_h=([0-9.]+)", phase)
    if m:
        return (1, float(m.group(1)), float(m.group(2)), phase)

    return (2, float("inf"), float("inf"), phase)

# ------------------------------------------------------------
# Discover phases automatically.
# Missing phases therefore do not cause any problem.
# ------------------------------------------------------------
phases = []
for entry in os.scandir(base_dir):
    if entry.is_dir():
        phase = entry.name
        if phase == "before_annealing" or phase == "bef_ann" or phase.startswith("annealing_T="):
            phases.append(phase)

phases = sorted(set(phases), key=phase_sort_key)

if not phases:
    raise SystemExit(
        f"ERROR: no annealing-phase directories found in:\n  {base_dir}"
    )

results = []

print()
print("============================================================")
print(" Total luminosity mean")
print("============================================================")
print(f" Sensor       : {sensor}")
print(f" Temperature  : T = {temperature_token} C")
print(f" Selected v_fin: {vfin_token} V")
print("============================================================")
print()

for phase in phases:

    # Expected structure:
    # DATA_irradiated_isolated_R=20/
    #   <phase>/
    #     A1_T=20_run=.../
    #       6luminosity/
    #         A1_T=20_<phase>_run=....csv

    phase_dir = os.path.join(base_dir, phase)

    pattern = os.path.join(
        phase_dir,
        f"{sensor}_T={temperature_token}_run=*",
        "6luminosity",
        f"{sensor}_T={temperature_token}_{phase}_run=*.csv",
    )

    files = sorted(glob.glob(pattern))

    # Small fallback for filenames/directories written with e.g. T=20.0.
    if not files:
        broad_pattern = os.path.join(
            phase_dir,
            f"{sensor}_T=*_run=*",
            "6luminosity",
            f"{sensor}_T=*_{phase}_run=*.csv",
        )

        candidate_files = sorted(glob.glob(broad_pattern))

        for path in candidate_files:
            run_dir = os.path.basename(
                os.path.dirname(os.path.dirname(path))
            )

            m = re.match(
                rf"^{re.escape(sensor)}_T=([-+0-9.eE]+)_run=",
                run_dir,
            )

            if not m:
                continue

            try:
                file_temperature = float(m.group(1))
            except ValueError:
                continue

            if math.isclose(
                file_temperature,
                temperature,
                rel_tol=0.0,
                abs_tol=1e-9,
            ):
                files.append(path)

    if not files:
        print(f"WARNING: {phase}: no matching CSV found -> phase skipped.")
        continue

    # Dictionary keyed by hotspot ID.
    # This guarantees that the same hotspot cannot silently be counted twice.
    hotspots = {}

    for path in files:
        with open(path, newline="") as f:
            reader = csv.DictReader(f)

            required = {"spot", "luminosity", "error", "v_fin"}
            missing = required.difference(reader.fieldnames or [])

            if missing:
                raise SystemExit(
                    "ERROR: missing required column(s) "
                    + ", ".join(sorted(missing))
                    + f"\nFile: {path}"
                )

            for line_number, row in enumerate(reader, start=2):

                try:
                    v_fin = float(row["v_fin"])
                except (TypeError, ValueError):
                    raise SystemExit(
                        f"ERROR: invalid v_fin at {path}:{line_number}"
                    )

                if not math.isclose(
                    v_fin,
                    target_vfin,
                    rel_tol=0.0,
                    abs_tol=1e-6,
                ):
                    continue

                spot_raw = row["spot"].strip()

                if spot_raw == "":
                    raise SystemExit(
                        f"ERROR: empty hotspot ID at {path}:{line_number}"
                    )

                # Normalize IDs such as "3" and "3.0" to the same key.
                try:
                    spot_value = float(spot_raw)
                    if spot_value.is_integer():
                        spot_id = str(int(spot_value))
                    else:
                        spot_id = f"{spot_value:g}"
                except ValueError:
                    spot_id = spot_raw

                try:
                    luminosity = float(row["luminosity"])
                    error = float(row["error"])
                except (TypeError, ValueError):
                    raise SystemExit(
                        f"ERROR: invalid luminosity/error at "
                        f"{path}:{line_number}"
                    )

                if not math.isfinite(luminosity):
                    raise SystemExit(
                        f"ERROR: non-finite luminosity at "
                        f"{path}:{line_number}"
                    )

                if not math.isfinite(error) or error < 0.0:
                    raise SystemExit(
                        f"ERROR: invalid statistical error at "
                        f"{path}:{line_number}"
                    )

                # Optional consistency check with the phase stored in the CSV.
                csv_phase = (row.get("phase") or "").strip()
                if csv_phase and csv_phase != phase:
                    raise SystemExit(
                        "ERROR: phase mismatch.\n"
                        f"Directory phase : {phase}\n"
                        f"CSV phase       : {csv_phase}\n"
                        f"File            : {path}\n"
                        f"Line            : {line_number}"
                    )

                if spot_id in hotspots:
                    previous = hotspots[spot_id]
                    raise SystemExit(
                        "ERROR: duplicate hotspot detected after the "
                        f"v_fin={vfin_token} selection.\n"
                        f"Phase      : {phase}\n"
                        f"Hotspot ID : {spot_id}\n"
                        f"First row  : {previous['file']}:{previous['line']}\n"
                        f"Second row : {path}:{line_number}\n"
                        "The phase was not processed, because counting the "
                        "same hotspot twice would bias the mean."
                    )

                hotspots[spot_id] = {
                    "luminosity": luminosity,
                    "error": error,
                    "file": path,
                    "line": line_number,
                }

    if not hotspots:
        print(
            f"WARNING: {phase}: no rows with v_fin={vfin_token} "
            "-> phase skipped."
        )
        continue

    n_spots = len(hotspots)

    total_luminosity = sum(
        item["luminosity"] for item in hotspots.values()
    )

    total_variance = sum(
        item["error"] ** 2 for item in hotspots.values()
    )

    mean_luminosity = total_luminosity / n_spots
    mean_error = math.sqrt(total_variance) / n_spots

    results.append(
        {
            "phase": phase,
            "mean_luminosity": mean_luminosity,
            "error": mean_error,
        }
    )

    print(
        f"{phase:28s}  "
        f"N = {n_spots:4d}  "
        f"<L> = {mean_luminosity:.10g} +/- {mean_error:.10g}"
    )

if not results:
    raise SystemExit(
        "\nERROR: no valid annealing phase produced a result."
    )

output_name = (
    f"{sensor}_T={temperature_token}_v={vfin_token}"
    "_total_luminosity_sum.csv"
)

output_path = os.path.join(output_dir, output_name)

with open(output_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["phase", "mean_luminosity", "error"])

    for row in results:
        writer.writerow(
            [
                row["phase"],
                f"{row['mean_luminosity']:.12g}",
                f"{row['error']:.12g}",
            ]
        )

print()
print("============================================================")
print("Output written to:")
print(f"  {output_path}")
print("============================================================")
PY
