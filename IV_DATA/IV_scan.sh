#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# IV_scan_makeiv.sh
#
# Standalone extraction pipeline reproducing the statistical logic of makeiv.C
# as closely as possible while producing one CSV per sensor.
#
# Location:
#   emmi/IV_DATA/IV_scan_makeiv.sh
#
# Input:
#   emmi/IV_DATA/A1/<phase>/*.ivscan.csv
#   emmi/IV_DATA/B1/<phase>/*.ivscan.csv
#
# Output:
#   emmi/IV_DATA/A1_IV_scan.csv
#   emmi/IV_DATA/B1_IV_scan.csv
#
# Output columns:
#   phase,temperature_C,voltage,current,current_error
#
# Important choices copied from makeiv.C:
#   - invert raw voltage sign: V = -V_raw
#   - invert raw current sign: I = -I_raw
#   - average repeated current measurements at fixed voltage
#   - TProfile-like default error:
#         sigma = sqrt(<I^2> - <I>^2)
#         error_on_mean = sigma / sqrt(N)
#   - subtract a zero-current level AFTER averaging
#   - combine IV and zero errors in quadrature
#   - DO NOT take fabs() after zero subtraction
#
# Zero level:
#   1) If a companion file containing "zero" and the same temperature token
#      (for example 290K) is found in the same phase directory, it is used
#      exactly like fnzero in makeiv.C: all third-column current measurements
#      are averaged into one zero-current level.
#   2) If no such file exists, the previous project convention is retained as
#      a fallback: mean currents at the 2nd, 3rd and 4th lowest distinct
#      voltages are averaged to estimate the zero-current level.
# ============================================================================

IV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 - "$IV_DIR" <<'PY'
import sys
import re
import math
import csv
from pathlib import Path
from collections import defaultdict

iv_dir = Path(sys.argv[1])
sensors = ("A1", "B1")


# ============================================================================
# Naming helpers
# ============================================================================

def get_phase(sensor, directory_name):
    name = directory_name
    low = name.lower()

    if "bef_ann" in low or "before_annealing" in low or "before" in low:
        return "bef_ann"

    m = re.search(
        r'(?:^|_)ann_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)(?:_|$)',
        name,
        re.IGNORECASE
    )
    if m:
        return f"annealing_T={m.group(1)}_h={m.group(2)}"

    m = re.search(
        r'annealing[_-]?T[=_-]?(\d+(?:\.\d+)?)[_-]?h[=_-]?(\d+(?:\.\d+)?)',
        name,
        re.IGNORECASE
    )
    if m:
        return f"annealing_T={m.group(1)}_h={m.group(2)}"

    phase = re.sub(rf'^{re.escape(sensor)}_?', '', name, flags=re.IGNORECASE)
    phase = re.sub(r'_IV_run=.*$', '', phase, flags=re.IGNORECASE)

    print(
        f"WARNING: cannot recognize phase from '{directory_name}'. "
        f"Using '{phase}'."
    )
    return phase


def temperature_token(filename):
    """Return the literal temperature token, e.g. 290K or 17C."""
    m = re.search(r'(\d+(?:\.\d+)?K)(?:_|\.|$)', filename, re.IGNORECASE)
    if m:
        return m.group(1)

    m = re.search(r'(\d+(?:\.\d+)?C)(?:_|\.|$)', filename, re.IGNORECASE)
    if m:
        return m.group(1)

    return None


def get_temperature_C(filename):
    m = re.search(r'(\d+(?:\.\d+)?)K(?:_|\.|$)', filename, re.IGNORECASE)
    if m:
        return int(round(float(m.group(1)) - 273.15))

    m = re.search(r'(\d+(?:\.\d+)?)C(?:_|\.|$)', filename, re.IGNORECASE)
    if m:
        return int(round(float(m.group(1))))

    return None


def phase_sort_key(phase):
    if phase in ("bef_ann", "before_annealing"):
        return (0, 0.0, 0.0, phase)

    m = re.match(r'annealing_T=([0-9.]+)_h=([0-9.]+)', phase)
    if m:
        return (1, float(m.group(1)), float(m.group(2)), phase)

    return (2, 0.0, 0.0, phase)


# ============================================================================
# makeiv/TProfile-like statistics
# ============================================================================

def profile_mean_and_error(values):
    """
    Reproduce the default unweighted TProfile statistics:

        mean = <x>
        sigma = sqrt(<x^2> - <x>^2)
        error = sigma / sqrt(N)

    This intentionally uses N (not N-1), because that is the convention used
    by TProfile::GetBinError() in the reference makeiv.C macro.
    """
    N = len(values)

    if N == 0:
        return None, None

    mean = sum(values) / N
    mean_sq = sum(x*x for x in values) / N

    variance = mean_sq - mean*mean

    # Floating-point roundoff can make a tiny negative value.
    if variance < 0.0 and abs(variance) < 1e-30:
        variance = 0.0

    variance = max(0.0, variance)

    sigma = math.sqrt(variance)
    error = sigma / math.sqrt(N)

    return mean, error


# ============================================================================
# Raw-file readers
# ============================================================================

def read_iv_file(filepath):
    """
    Equivalent sign convention to:
        graphutils::invertX(givscan)
        graphutils::invertY(givscan)

    Returns:
        voltage -> [current measurements]

    Voltage is quantized to 1 mV to reproduce the binning of:
        TProfile(..., 100001, -0.0005, 100.0005)
    whose bin width is 0.001 V.
    """
    data = defaultdict(list)

    with filepath.open("r") as fin:
        for line in fin:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 3:
                continue

            try:
                raw_voltage = float(fields[1])
                raw_current = float(fields[2])
            except ValueError:
                continue

            voltage = -raw_voltage
            current = -raw_current

            # TProfile binning has 1 mV bins centered at integer millivolts.
            voltage = round(voltage, 3)

            if not math.isfinite(voltage) or not math.isfinite(current):
                continue

            data[voltage].append(current)

    return data


def read_zero_file(filepath):
    """
    Match makeiv.C zero reader:
        TGraph(fnzero, "%lg %*lg %lg")
        invertY

    The second column is intentionally ignored.
    """
    currents = []

    with filepath.open("r") as fin:
        for line in fin:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 3:
                continue

            try:
                raw_current = float(fields[2])
            except ValueError:
                continue

            current = -raw_current

            if math.isfinite(current):
                currents.append(current)

    return currents


def find_zero_file(phase_dir, iv_file):
    """
    Look for a companion zero-level file in the same phase directory.
    It must contain 'zero' in its name and, when available, the same
    temperature token as the IV file.
    """
    token = temperature_token(iv_file.name)

    candidates = []

    for p in sorted(phase_dir.iterdir()):
        if not p.is_file():
            continue

        low = p.name.lower()

        if "zero" not in low:
            continue

        if p.suffix.lower() == ".root":
            continue

        if token is not None and token.lower() not in low:
            continue

        candidates.append(p)

    if not candidates:
        return None

    if len(candidates) > 1:
        print(
            f"WARNING: multiple zero files for {iv_file.name}; "
            f"using {candidates[0].name}"
        )

    return candidates[0]


# ============================================================================
# Zero-current estimate
# ============================================================================

def fallback_zero_from_low_voltage(points_by_voltage):
    """
    Fallback retained from the previous pipeline:
      - take 2nd, 3rd, 4th lowest distinct voltages
      - calculate mean current at each fixed voltage
      - average those three means
      - propagate their mean errors
    """
    voltages = sorted(points_by_voltage.keys())

    if len(voltages) < 4:
        return None

    selected = voltages[1:4]
    components = []

    for voltage in selected:
        mean_I, error_I = profile_mean_and_error(points_by_voltage[voltage])

        if mean_I is None:
            return None

        components.append((voltage, mean_I, error_I))

    zero = sum(x[1] for x in components) / 3.0
    ezero = math.sqrt(sum(x[2]**2 for x in components)) / 3.0

    return zero, ezero, components


# ============================================================================
# Process sensors
# ============================================================================

for sensor in sensors:

    sensor_dir = iv_dir / sensor
    output_file = iv_dir / f"{sensor}_IV_scan.csv"

    print()
    print("=" * 72)
    print(f"SENSOR: {sensor}")
    print(f"Input : {sensor_dir}")
    print(f"Output: {output_file}")
    print("=" * 72)

    if not sensor_dir.is_dir():
        print(f"WARNING: missing sensor directory: {sensor_dir}")
        continue

    output_rows = []
    n_iv_files = 0
    n_zero_external = 0
    n_zero_fallback = 0

    for phase_dir in sorted(sensor_dir.iterdir()):

        if not phase_dir.is_dir():
            continue

        phase = get_phase(sensor, phase_dir.name)

        iv_files = [
            p for p in sorted(phase_dir.glob("*.ivscan.csv"))
            if "zoom" not in p.name.lower()
        ]

        for iv_file in iv_files:

            n_iv_files += 1

            temperature_C = get_temperature_C(iv_file.name)

            if temperature_C is None:
                print(
                    f"WARNING: cannot extract temperature from "
                    f"{iv_file.name}; skipping."
                )
                continue

            points_by_voltage = read_iv_file(iv_file)

            if not points_by_voltage:
                print(f"WARNING: no valid points in {iv_file.name}; skipping.")
                continue

            # ----------------------------------------------------------------
            # Zero level: prefer an external zero file, as makeiv.C does.
            # ----------------------------------------------------------------

            zero_file = find_zero_file(phase_dir, iv_file)

            if zero_file is not None:
                zero_values = read_zero_file(zero_file)
                zero, ezero = profile_mean_and_error(zero_values)

                if zero is None:
                    print(
                        f"WARNING: zero file {zero_file.name} has no valid "
                        f"measurements; using low-voltage fallback."
                    )
                    zero_file = None
                else:
                    n_zero_external += 1
                    zero_source = f"external:{zero_file.name}"

            if zero_file is None:
                fallback = fallback_zero_from_low_voltage(points_by_voltage)

                if fallback is None:
                    print(
                        f"ERROR: cannot estimate zero for {iv_file.name}: "
                        f"fewer than four distinct voltages."
                    )
                    continue

                zero, ezero, components = fallback
                n_zero_fallback += 1
                zero_source = "fallback:2nd-4th-lowest-voltage"

            print()
            print("-" * 72)
            print(f"{sensor} | {phase} | T={temperature_C} C")
            print(f"file        : {iv_file.name}")
            print(f"zero source : {zero_source}")
            print(f"zero        : {zero:.12e} +/- {ezero:.12e} A")

            # ----------------------------------------------------------------
            # Build makeiv-like curve.
            # ----------------------------------------------------------------

            for voltage in sorted(points_by_voltage.keys()):

                mean_I, profile_error = profile_mean_and_error(
                    points_by_voltage[voltage]
                )

                if mean_I is None:
                    continue

                # makeiv.C explicitly skips profile bins with zero error.
                if profile_error == 0.0:
                    continue

                corrected_current = mean_I - zero

                corrected_error = math.sqrt(
                    profile_error*profile_error + ezero*ezero
                )

                output_rows.append(
                    (
                        phase,
                        temperature_C,
                        voltage,
                        corrected_current,
                        corrected_error
                    )
                )

    output_rows.sort(
        key=lambda row: (
            phase_sort_key(row[0]),
            row[1],
            row[2]
        )
    )

    with output_file.open("w", newline="") as fout:
        writer = csv.writer(fout)

        writer.writerow([
            "phase",
            "temperature_C",
            "voltage",
            "current",
            "current_error"
        ])

        for phase, temperature_C, voltage, current, current_error in output_rows:
            writer.writerow([
                phase,
                temperature_C,
                f"{voltage:.6f}",
                f"{current:.12e}",
                f"{current_error:.12e}"
            ])

    print()
    print(f"Completed {sensor}")
    print(f"  IV files processed       : {n_iv_files}")
    print(f"  external zero files used : {n_zero_external}")
    print(f"  low-voltage fallbacks    : {n_zero_fallback}")
    print(f"  output rows              : {len(output_rows)}")
    print(f"  output CSV               : {output_file}")

print()
print("=" * 72)
print("makeiv-like IV extraction completed")
print("=" * 72)
PY
