#!/usr/bin/env bash
set -euo pipefail

# Usage: ./IV_I_vs_T.sh A1|B1
# Script location: emmi/IV_DATA/
# Output: emmi/IV_DATA/<SENSOR>_I_vs_T.csv

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 SENSOR"
    echo "  SENSOR = A1 or B1"
    exit 1
fi

SENSOR="$1"
case "$SENSOR" in
    A1|B1) ;;
    *) echo "ERROR: sensor must be A1 or B1."; exit 1 ;;
esac

IV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SENSOR_DIR="${IV_DIR}/${SENSOR}"
OUTFILE="${IV_DIR}/${SENSOR}_I_vs_T.csv"

if [[ ! -d "$SENSOR_DIR" ]]; then
    echo "ERROR: directory not found: $SENSOR_DIR"
    exit 1
fi

echo "============================================================"
echo " IV current vs temperature analysis"
echo " Sensor:   $SENSOR"
echo " Input:    $SENSOR_DIR"
echo " Output:   $OUTFILE"
echo "============================================================"

python3 - "$SENSOR" "$SENSOR_DIR" "$OUTFILE" <<'PY'
import sys, re, math, csv, statistics
from pathlib import Path
from collections import defaultdict

sensor = sys.argv[1]
sensor_dir = Path(sys.argv[2])
output_file = Path(sys.argv[3])

TARGET_VOLTAGES = {
    "A1": {17: 56.12, 19: 56.24, 21: 56.35, 23: 56.46},
    "B1": {17: 53.75, 19: 53.85, 21: 53.96, 23: 54.06},
}
targets = TARGET_VOLTAGES[sensor]


def get_phase(directory_name):
    name = directory_name
    lower = name.lower()
    if "before" in lower or "bef_ann" in lower:
        return "bef_ann"

    m = re.search(r'(?:^|_)ann_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)(?:_|$)', name, re.I)
    if m:
        return f"annealing_T={m.group(1)}_h={m.group(2)}"

    m = re.search(r'annealing[_-]?T[=_-]?(\d+(?:\.\d+)?)[_-]?h[=_-]?(\d+(?:\.\d+)?)', name, re.I)
    if m:
        return f"annealing_T={m.group(1)}_h={m.group(2)}"

    phase = re.sub(rf'^{re.escape(sensor)}_?', '', name, flags=re.I)
    phase = re.sub(r'_IV_run=.*$', '', phase, flags=re.I)
    print(f"WARNING: could not recognize phase from '{directory_name}'. Using '{phase}'.")
    return phase


def get_temperature_C(filename):
    m = re.search(r'(\d+(?:\.\d+)?)K(?:_|\.|$)', filename, re.I)
    if m:
        return int(round(float(m.group(1)) - 273.15))

    m = re.search(r'(\d+(?:\.\d+)?)C(?:_|\.|$)', filename, re.I)
    if m:
        return int(round(float(m.group(1))))

    return None


def phase_sort_key(phase):
    if phase in ("bef_ann", "before_annealing"):
        return (0, 0.0, 0.0)
    m = re.match(r'annealing_T=([0-9.]+)_h=([0-9.]+)', phase)
    if m:
        return (1, float(m.group(1)), float(m.group(2)))
    return (2, 0.0, 0.0)


def mean_and_error(values):
    n = len(values)
    if n == 0:
        return None, None
    mean = statistics.fmean(values)
    if n > 1:
        err = statistics.stdev(values) / math.sqrt(n)
    else:
        err = float("nan")
    return mean, err


def build_voltage_points(voltage_data):
    points = []
    for voltage in sorted(voltage_data):
        values = voltage_data[voltage]
        mean_I, error_I = mean_and_error(values)
        if mean_I is not None:
            points.append((voltage, mean_I, error_I, len(values)))
    return points


def interpolate_current(voltage_data, target_voltage):
    points = build_voltage_points(voltage_data)
    if not points:
        return None

    # Exact match, if ever present.
    for voltage, mean_I, error_I, n in points:
        if abs(voltage - target_voltage) <= 1e-9:
            return {
                "current": mean_I, "error": error_I,
                "v0": voltage, "v1": voltage,
                "I0": mean_I, "I1": mean_I,
                "e0": error_I, "e1": error_I,
                "N0": n, "N1": n, "fraction": 0.0,
            }

    lower = None
    upper = None
    for point in points:
        v = point[0]
        if v < target_voltage:
            lower = point
        elif v > target_voltage:
            upper = point
            break

    if lower is None or upper is None:
        return None

    v0, I0, e0, N0 = lower
    v1, I1, e1, N1 = upper
    if abs(v1 - v0) < 1e-15:
        return None

    t = (target_voltage - v0) / (v1 - v0)
    I = I0 + t * (I1 - I0)

    if math.isfinite(e0) and math.isfinite(e1):
        err = math.sqrt(((1.0 - t) * e0)**2 + (t * e1)**2)
    else:
        err = float("nan")

    return {
        "current": I, "error": err,
        "v0": v0, "v1": v1,
        "I0": I0, "I1": I1,
        "e0": e0, "e1": e1,
        "N0": N0, "N1": N1, "fraction": t,
    }


def calculate_offset(voltage_data):
    """
    Offset = mean of the mean currents at the 2nd, 3rd and 4th
    lowest DISTINCT measured voltage values.
    """
    points = build_voltage_points(voltage_data)
    if len(points) < 4:
        return None

    # index 0 = lowest; indices 1,2,3 = 2nd,3rd,4th lowest
    selected = points[1:4]
    offset = sum(p[1] for p in selected) / 3.0

    errors = [p[2] for p in selected]
    if all(math.isfinite(e) for e in errors):
        offset_error = math.sqrt(sum(e * e for e in errors)) / 3.0
    else:
        offset_error = float("nan")

    return {
        "offset": offset,
        "error": offset_error,
        "points": selected,
    }


# scans[(phase, temperature_C)][voltage] = [current samples]
scans = defaultdict(lambda: defaultdict(list))
n_files_found = 0
n_files_used = 0

for phase_dir in sorted(sensor_dir.iterdir()):
    if not phase_dir.is_dir():
        continue

    phase = get_phase(phase_dir.name)

    for filepath in sorted(phase_dir.glob("*.ivscan.csv")):
        n_files_found += 1

        if "zoom" in filepath.name.lower():
            print(f"Skipping zoom file: {filepath.name}")
            continue

        temperature_C = get_temperature_C(filepath.name)
        if temperature_C is None:
            print(f"WARNING: cannot determine temperature from {filepath.name}. Skipping.")
            continue

        if temperature_C not in targets:
            print(
                f"WARNING: no target voltage defined for {sensor}, "
                f"T = {temperature_C} C. Skipping {filepath.name}."
            )
            continue

        print()
        print(f"Reading {filepath.name}")
        print(f"  phase = {phase}")
        print(f"  T     = {temperature_C} C")
        print(f"  target voltage = {targets[temperature_C]:.2f} V")

        n_rows = 0
        with filepath.open("r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                columns = line.split()
                if len(columns) < 3:
                    continue

                try:
                    voltage = abs(float(columns[1]))
                    current = abs(float(columns[2]))
                except ValueError:
                    continue

                # Group repeated samples at the same voltage setting.
                voltage = round(voltage, 6)
                scans[(phase, temperature_C)][voltage].append(current)
                n_rows += 1

        if n_rows > 0:
            n_files_used += 1
        else:
            print("WARNING: no valid IV rows found.")


results = []

for (phase, temperature_C), voltage_data in scans.items():
    target_voltage = targets[temperature_C]

    interpolation = interpolate_current(voltage_data, target_voltage)
    offset_result = calculate_offset(voltage_data)

    print()
    print("------------------------------------------------------------")
    print(f"{sensor} | {phase} | T = {temperature_C} C")
    print(f"Target V = {target_voltage:.2f} V")

    if interpolation is None:
        available = sorted(voltage_data)
        print("ERROR: cannot interpolate target voltage.")
        if available:
            print(f"Available voltage range: {available[0]:.3f} - {available[-1]:.3f} V")
        continue

    if offset_result is None:
        print("ERROR: cannot calculate current offset.")
        print("At least 4 distinct voltage values are required.")
        continue

    print()
    print("Interpolation points:")
    print(
        f"  V0 = {interpolation['v0']:.3f} V"
        f"   I0 = {interpolation['I0']:.6e} A"
        f"   err0 = {interpolation['e0']:.6e} A"
        f"   N0 = {interpolation['N0']}"
    )
    print(
        f"  V1 = {interpolation['v1']:.3f} V"
        f"   I1 = {interpolation['I1']:.6e} A"
        f"   err1 = {interpolation['e1']:.6e} A"
        f"   N1 = {interpolation['N1']}"
    )
    print(f"Interpolation fraction t = {interpolation['fraction']:.6f}")
    print(
        f"Interpolated current before offset = "
        f"{interpolation['current']:.6e} +/- {interpolation['error']:.6e} A"
    )

    print()
    print("Offset calculation:")
    ordinal = ["2nd", "3rd", "4th"]
    for label, (voltage, mean_I, error_I, n) in zip(ordinal, offset_result["points"]):
        print(
            f"  {label} lowest V = {voltage:.3f} V"
            f"   mean I = {mean_I:.6e} A"
            f"   err = {error_I:.6e} A"
            f"   N = {n}"
        )

    offset_current = offset_result["offset"]
    offset_error = offset_result["error"]
    print(f"Offset = {offset_current:.6e} +/- {offset_error:.6e} A")

    corrected_current = interpolation["current"] - offset_current

    if math.isfinite(interpolation["error"]) and math.isfinite(offset_error):
        corrected_error = math.sqrt(interpolation["error"]**2 + offset_error**2)
    else:
        corrected_error = float("nan")

    print()
    print(f"Corrected current = {corrected_current:.6e} +/- {corrected_error:.6e} A")

    results.append((phase, temperature_C, corrected_current, corrected_error))


if not results:
    print()
    print("ERROR: no corrected current values produced.")
    sys.exit(1)

results.sort(key=lambda x: (phase_sort_key(x[0]), x[1]))

with output_file.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["phase", "temperature_C", "current", "current_error"])

    for phase, temperature_C, current, current_error in results:
        writer.writerow([
            phase,
            temperature_C,
            f"{current:.12e}",
            f"{current_error:.12e}",
        ])

print()
print("============================================================")
print(" Analysis completed")
print(f" IV files found:     {n_files_found}")
print(f" IV files used:      {n_files_used}")
print(f" Output rows:        {len(results)}")
print(f" Output file:        {output_file}")
print("============================================================")

PY
