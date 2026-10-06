#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./hotspot_photos.sh SENSOR GLOBAL_ID OPERATION

SENSOR:
  A1 | B1

OPERATION:
  annealing_scan | temperature_scan | voltage_scan

Examples:
  ./hotspot_photos.sh A1 0 annealing_scan
  ./hotspot_photos.sh A1 0 temperature_scan
  ./hotspot_photos.sh B1 12 voltage_scan
EOF
}

if [ "$#" -ne 3 ]; then
    usage
    exit 1
fi

SENSOR="$(printf '%s' "$1" | tr '[:lower:]' '[:upper:]')"
GLOBAL_ID="$2"
OPERATION="$3"

case "$SENSOR" in
    A1)
        FIXED_VOVER="5"
        ;;
    B1)
        FIXED_VOVER="3"
        ;;
    *)
        echo "ERROR: SENSOR must be A1 or B1." >&2
        exit 1
        ;;
esac

case "$GLOBAL_ID" in
    ''|*[!0-9]*)
        echo "ERROR: GLOBAL_ID must be a non-negative integer." >&2
        exit 1
        ;;
esac

case "$OPERATION" in
    annealing_scan|temperature_scan|voltage_scan)
        ;;
    *)
        echo "ERROR: OPERATION must be annealing_scan, temperature_scan, or voltage_scan." >&2
        exit 1
        ;;
esac

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EMMI_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_ROOT="$EMMI_DIR/DATA_irradiated_isolated_R=20"
PYTHON_BIN="${PYTHON_BIN:-python}"

# Accept either of the two sensible locations for find_centers_save.py.
if [ -f "$EMMI_DIR/find_centers/find_centers_save.py" ]; then
    FIND_CENTERS="$EMMI_DIR/find_centers/find_centers_save.py"
elif [ -f "$SCRIPT_DIR/find_centers_save.py" ]; then
    FIND_CENTERS="$SCRIPT_DIR/find_centers_save.py"
else
    echo "ERROR: find_centers_save.py not found." >&2
    echo "Checked:" >&2
    echo "  $EMMI_DIR/find_centers/find_centers_save.py" >&2
    echo "  $SCRIPT_DIR/find_centers_save.py" >&2
    exit 1
fi

case "$OPERATION" in
    annealing_scan|voltage_scan)
        CONDITION_DIR="${SENSOR}_T=20"
        MAPPING_FILE="$DATA_ROOT/merged_files/$CONDITION_DIR/${SENSOR}_T_global_mapping.csv"
        GLOBAL_FILE="$DATA_ROOT/merged_files/$CONDITION_DIR/${SENSOR}_T_all_phases_global_ID.csv"
        ;;
    temperature_scan)
        CONDITION_DIR="${SENSOR}_v=${FIXED_VOVER}"
        MAPPING_FILE="$DATA_ROOT/merged_files/$CONDITION_DIR/${SENSOR}_v_global_mapping.csv"
        GLOBAL_FILE="$DATA_ROOT/merged_files/$CONDITION_DIR/${SENSOR}_v_all_phases_global_ID.csv"
        ;;
esac

for required_file in "$MAPPING_FILE" "$GLOBAL_FILE" "$FIND_CENTERS"; do
    if [ ! -f "$required_file" ]; then
        echo "ERROR: required file not found: $required_file" >&2
        exit 1
    fi
done

OUTPUT_DIR="$SCRIPT_DIR/$SENSOR/$OPERATION/spot_${GLOBAL_ID}"
mkdir -p "$OUTPUT_DIR"

printf '%s\n' "============================================================"
printf 'Hotspot photo pipeline\n'
printf 'Sensor        : %s\n' "$SENSOR"
printf 'Global ID     : %s\n' "$GLOBAL_ID"
printf 'Operation     : %s\n' "$OPERATION"
printf 'Mapping CSV   : %s\n' "$MAPPING_FILE"
printf 'Global CSV    : %s\n' "$GLOBAL_FILE"
printf 'Finder        : %s\n' "$FIND_CENTERS"
printf 'Output dir    : %s\n' "$OUTPUT_DIR"
printf '%s\n' "============================================================"

# Keep the path/CSV orchestration in Python. This avoids Bash parsing problems
# with phases, run numbers and floating-point values while the actual source
# finding remains entirely inside find_centers_save.py.
"$PYTHON_BIN" - \
    "$DATA_ROOT" \
    "$MAPPING_FILE" \
    "$GLOBAL_FILE" \
    "$FIND_CENTERS" \
    "$OUTPUT_DIR" \
    "$SENSOR" \
    "$GLOBAL_ID" \
    "$OPERATION" \
    "$FIXED_VOVER" <<'PY'
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import sys

try:
    import tifffile
except ImportError as exc:
    raise SystemExit(
        "ERROR: tifffile is required in the active Python environment"
    ) from exc

(
    data_root_s,
    mapping_file_s,
    global_file_s,
    find_centers_s,
    output_dir_s,
    sensor,
    global_id_s,
    operation,
    fixed_vover_s,
) = sys.argv[1:]

data_root = Path(data_root_s)
mapping_file = Path(mapping_file_s)
global_file = Path(global_file_s)
find_centers = Path(find_centers_s)
output_dir = Path(output_dir_s)
global_id = int(global_id_s)
fixed_vover = float(fixed_vover_s)


def as_float(value, field, context=""):
    try:
        result = float(value)
    except (TypeError, ValueError):
        raise RuntimeError(
            f"invalid numeric value for {field}: {value!r} {context}"
        )
    if not math.isfinite(result):
        raise RuntimeError(
            f"non-finite value for {field}: {value!r} {context}"
        )
    return result


def as_int(value, field, context=""):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        raise RuntimeError(
            f"invalid integer value for {field}: {value!r} {context}"
        )


def num_equal(a, b, tol=1.0e-6):
    return abs(float(a) - float(b)) < tol


def fmt_num(value):
    value = float(value)
    if abs(value - round(value)) < 1.0e-9:
        return str(int(round(value)))
    return f"{value:g}"


def normalize_phase(phase):
    phase = phase.strip()
    if phase == "before_annealing":
        return phase

    match = re.fullmatch(r"annealing_T=(\d+)_h=(\d+)", phase)
    if match:
        return f"annealing_{match.group(1)}_{match.group(2)}"

    match = re.fullmatch(r"annealing_(\d+)_(\d+)", phase)
    if match:
        return phase

    return phase


def old_phase_name(phase):
    normalized = normalize_phase(phase)
    match = re.fullmatch(r"annealing_(\d+)_(\d+)", normalized)
    if match:
        return f"annealing_T={match.group(1)}_h={match.group(2)}"
    return normalized


def pretty_phase(phase):
    normalized = normalize_phase(phase)
    if normalized == "before_annealing":
        return "Before annealing"
    match = re.fullmatch(r"annealing_(\d+)_(\d+)", normalized)
    if match:
        return f"Annealing {match.group(1)} °C, {match.group(2)} h"
    return phase


def phase_sort_key(phase):
    normalized = normalize_phase(phase)
    if normalized == "before_annealing":
        return (0, 0, 0)
    match = re.fullmatch(r"annealing_(\d+)_(\d+)", normalized)
    if match:
        return (1, int(match.group(1)), int(match.group(2)))
    return (2, normalized, 0)


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = set(reader.fieldnames or [])
    return rows, fields


global_rows, global_fields = read_csv(global_file)
mapping_rows, mapping_fields = read_csv(mapping_file)

required_global = {
    "spot",
    "x",
    "y",
    "T",
    "v",
    "v_fin",
    "phase",
    "detection",
    "run_number",
}
missing = required_global - global_fields
if missing:
    raise SystemExit(
        f"ERROR: missing columns in {global_file}: {', '.join(sorted(missing))}"
    )

required_mapping = {"phase", "constant", "run_number", "local_spot", "spot"}
missing = required_mapping - mapping_fields
if missing:
    raise SystemExit(
        f"ERROR: missing columns in {mapping_file}: {', '.join(sorted(missing))}"
    )

# ------------------------------------------------------------------
# Select the rows of *_all_phases_global_ID.csv belonging to this scan.
# ------------------------------------------------------------------
spot_rows = []
for row in global_rows:
    try:
        spot = int(float(row["spot"]))
    except (TypeError, ValueError):
        continue
    if spot == global_id:
        spot_rows.append(row)

if not spot_rows:
    raise SystemExit(
        f"ERROR: global hotspot {global_id} not found in {global_file}"
    )

selected = []
if operation == "annealing_scan":
    for row in spot_rows:
        v_fin = as_float(row["v_fin"], "v_fin")
        if num_equal(v_fin, fixed_vover):
            selected.append(row)
    selected.sort(key=lambda row: phase_sort_key(row["phase"]))

elif operation == "voltage_scan":
    selected = [
        row
        for row in spot_rows
        if normalize_phase(row["phase"]) == "before_annealing"
    ]
    selected.sort(key=lambda row: as_float(row["v_fin"], "v_fin"))

elif operation == "temperature_scan":
    selected = [
        row
        for row in spot_rows
        if normalize_phase(row["phase"]) == "before_annealing"
    ]
    selected.sort(key=lambda row: as_float(row["T"], "T"))

if not selected:
    raise SystemExit(
        f"ERROR: no rows selected for {sensor} hotspot {global_id} "
        f"operation {operation}"
    )

# Check uniqueness of the intended scan points.
seen_keys = set()
for row in selected:
    if operation == "annealing_scan":
        key = normalize_phase(row["phase"])
    elif operation == "voltage_scan":
        key = as_float(row["v_fin"], "v_fin")
    else:
        key = as_float(row["T"], "T")
    if key in seen_keys:
        raise SystemExit(
            f"ERROR: duplicate scan point {key!r} in {global_file} for hotspot {global_id}"
        )
    seen_keys.add(key)

# The mapping's `constant` is the fixed dataset condition:
#   *_T=20  -> constant = 20
#   *_v=5/3 -> constant = 5/3
mapping_constant = 20.0 if operation in {"annealing_scan", "voltage_scan"} else fixed_vover


def lookup_local_spot(phase, run_number):
    candidates = []
    for row in mapping_rows:
        try:
            spot = int(float(row["spot"]))
        except (TypeError, ValueError):
            continue

        if spot != global_id:
            continue
        if normalize_phase(row["phase"]) != normalize_phase(phase):
            continue

        constant_raw = (row.get("constant") or "").strip()
        constant_ok = False
        try:
            constant_ok = num_equal(float(constant_raw), mapping_constant)
        except ValueError:
            # Backward compatibility with mappings that store the fixed
            # condition as a label rather than as its numeric value.
            if operation in {"annealing_scan", "voltage_scan"}:
                constant_ok = constant_raw.lower() == "t"
            else:
                constant_ok = constant_raw.lower() == "v"

        if not constant_ok:
            continue

        # Prefer the exact run when the mapping contains it.
        map_run = (row.get("run_number") or "").strip()
        if map_run and run_number and map_run != run_number:
            continue

        candidates.append(as_int(row["local_spot"], "local_spot"))

    unique = sorted(set(candidates))
    if len(unique) == 1:
        return unique[0]
    if len(unique) == 0:
        return None
    raise RuntimeError(
        f"ambiguous local_spot for global hotspot {global_id}, phase {phase}: {unique}"
    )


def resolve_phase_dir(phase):
    candidates = [phase, normalize_phase(phase), old_phase_name(phase)]
    seen = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        path = data_root / candidate
        if path.is_dir():
            return path
    raise RuntimeError(
        f"data directory not found for phase {phase!r}; checked "
        + ", ".join(str(data_root / c) for c in seen)
    )


def parse_tif_conditions(path):
    match = re.search(
        r"_T=([-+0-9.eE]+)_v=([-+0-9.eE]+)_data=diff_processed\.(?:tif|tiff)$",
        path.name,
        re.IGNORECASE,
    )
    if not match:
        return None
    return float(match.group(1)), float(match.group(2))


def find_tif(row):
    phase = row["phase"].strip()
    run_number = row["run_number"].strip()
    temp = as_float(row["T"], "T", f"in phase {phase}")
    voltage = as_float(row["v"], "v", f"in phase {phase}")

    phase_dir = resolve_phase_dir(phase)

    if operation == "temperature_scan":
        dataset_dir = f"{sensor}_v={fmt_num(fixed_vover)}_run={run_number}"
    else:
        dataset_dir = f"{sensor}_T=20_run={run_number}"

    processed_dir = phase_dir / dataset_dir / "2processed"
    if not processed_dir.is_dir():
        raise RuntimeError(f"processed directory not found: {processed_dir}")

    candidates = []
    for pattern in ("*.tif", "*.tiff", "*.TIF", "*.TIFF"):
        for path in processed_dir.glob(pattern):
            # EXACTLY diff_processed; never diffe_processed.
            if not re.search(r"_data=diff_processed\.(?:tif|tiff)$", path.name, re.I):
                continue
            conditions = parse_tif_conditions(path)
            if conditions is None:
                continue
            tif_temp, tif_voltage = conditions
            if num_equal(tif_temp, temp) and num_equal(tif_voltage, voltage):
                candidates.append(path)

    unique = sorted(set(candidates))
    if len(unique) == 0:
        raise RuntimeError(
            f"no diff_processed TIF at T={temp:g}, v={voltage:g} in {processed_dir}"
        )
    if len(unique) > 1:
        raise RuntimeError(
            "more than one matching diff_processed TIF:\n  "
            + "\n  ".join(str(path) for path in unique)
        )
    return unique[0]


def root_to_python(path, x_root, y_root):
    image = tifffile.imread(path)
    if image.ndim != 2:
        raise RuntimeError(f"expected 2D TIF, got {image.shape} for {path}")
    ny, nx = image.shape
    x_python = float(x_root)
    y_python = float(ny) - float(y_root)
    if not (0 <= x_python < nx and 0 <= y_python < ny):
        raise RuntimeError(
            f"converted focus_area x={x_python:g}, y={y_python:g} is outside "
            f"image {nx}x{ny}; ROOT input was x={x_root}, y={y_root}"
        )
    return x_python, y_python, nx, ny


def output_name(row):
    if operation == "annealing_scan":
        return normalize_phase(row["phase"])
    if operation == "voltage_scan":
        return f"v_over_{fmt_num(as_float(row['v_fin'], 'v_fin'))}"
    return f"T_{fmt_num(as_float(row['T'], 'T'))}C"


def build_title(row):
    temp = as_float(row["T"], "T")
    v_fin = as_float(row["v_fin"], "v_fin")
    return (
        f"{sensor} — global hotspot {global_id} — {pretty_phase(row['phase'])}\n"
        f"T = {fmt_num(temp)} °C, v_over = {fmt_num(v_fin)} V"
    )

# ------------------------------------------------------------------
# Build one explicit entry per image.
# ------------------------------------------------------------------
entries = []
last_detected_root_xy = None

for row in selected:
    phase = row["phase"].strip()
    run_number = row["run_number"].strip()
    detection = (row.get("detection") or "").strip().lower()
    x_root = as_float(row["x"], "x", f"in phase {phase}")
    y_root = as_float(row["y"], "y", f"in phase {phase}")
    tif = find_tif(row)

    if detection == "present":
        local_spot = lookup_local_spot(phase, run_number)
        if local_spot is None:
            raise RuntimeError(
                f"hotspot {global_id} is marked present in phase {phase}, but no "
                f"local_spot was found in {mapping_file}"
            )
        mode = "focus"
        focus_args = (local_spot,)
        last_detected_root_xy = (x_root, y_root)
    else:
        # For an absent hotspot, use the coordinates of the last detection.
        # The all-phases CSV retains coordinates also for absent rows, so if
        # there has not yet been a previous detection we use those stored
        # coordinates rather than aborting the scan.
        if last_detected_root_xy is None:
            last_detected_root_xy = (x_root, y_root)
            print(
                f"WARNING: no earlier detected point for {phase}; using the "
                f"stored global coordinates x={x_root:g}, y={y_root:g}",
                file=sys.stderr,
            )

        x_last_root, y_last_root = last_detected_root_xy
        x_python, y_python, nx, ny = root_to_python(
            tif, x_last_root, y_last_root
        )
        mode = "focus_area"
        focus_args = (x_python, y_python)

    entries.append(
        {
            "row": row,
            "phase": phase,
            "run": run_number,
            "detection": detection,
            "tif": tif,
            "mode": mode,
            "focus_args": focus_args,
            "output": output_dir / output_name(row),
            "title": build_title(row),
        }
    )

# ------------------------------------------------------------------
# Helper to invoke find_centers_save.py.
# ------------------------------------------------------------------
env = os.environ.copy()
env["MPLBACKEND"] = "Agg"


def make_command(entry, measure=False, vmin=None, vmax=None):
    cmd = [
        sys.executable,
        str(find_centers),
        "--input",
        str(entry["tif"]),
    ]

    if entry["mode"] == "focus":
        cmd += ["--focus", str(entry["focus_args"][0])]
    else:
        cmd += [
            "--focus_area",
            f"{entry['focus_args'][0]:.12g}",
            f"{entry['focus_args'][1]:.12g}",
        ]

    if measure:
        cmd.append("--measure_max")
    else:
        cmd += [
            "--output_focus",
            str(entry["output"]),
            "--focus_title",
            entry["title"],
            "--vmin",
            f"{vmin:.12g}",
            "--vmax",
            f"{vmax:.12g}",
        ]

    return cmd

# ------------------------------------------------------------------
# First pass: calculate ONE common display vmax for this hotspot + scan.
# ------------------------------------------------------------------
maxima = []
print("\n--- First pass: common colour scale ---")
for entry in entries:
    print(
        f"  {entry['phase']}: detection={entry['detection'] or 'unknown'}, "
        f"mode={entry['mode']}, tif={entry['tif'].name}"
    )
    result = subprocess.run(
        make_command(entry, measure=True),
        env=env,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        raise RuntimeError(
            f"find_centers_save.py failed while measuring scale for {entry['phase']}"
        )

    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if not lines:
        raise RuntimeError(
            f"no scale value returned by find_centers_save.py for {entry['phase']}"
        )
    try:
        maxima.append(float(lines[-1]))
    except ValueError:
        raise RuntimeError(
            f"invalid scale output for {entry['phase']}: {lines[-1]!r}"
        )

vmin = 0.0
vmax = max(maxima)
if not math.isfinite(vmax):
    raise RuntimeError("common colour-scale vmax is not finite")
if vmax <= vmin:
    vmax = 1.0

scale_file = output_dir / "display_scale.csv"
with scale_file.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["vmin", "vmax", "n_images"])
    writer.writerow([f"{vmin:.12g}", f"{vmax:.12g}", len(entries)])

print(f"Common display scale: vmin={vmin:g}, vmax={vmax:g}")
print(f"Scale saved to: {scale_file}")

# ------------------------------------------------------------------
# Second pass: save every image with the SAME scale.
# ------------------------------------------------------------------
print("\n--- Second pass: save focused images ---")
for index, entry in enumerate(entries, start=1):
    row = entry["row"]
    print("------------------------------------------------------------")
    print(f"[{index}/{len(entries)}] Phase       : {entry['phase']}")
    print(f"Detection            : {entry['detection'] or 'unknown'}")
    print(f"Run                  : {entry['run']}")
    print(f"TIF                  : {entry['tif']}")
    if entry["mode"] == "focus":
        print(f"Focus mode           : --focus local_spot={entry['focus_args'][0]}")
    else:
        print(
            "Focus mode           : --focus_area "
            f"x={entry['focus_args'][0]:.3f}, y={entry['focus_args'][1]:.3f} "
            "[Python coordinates from last detection]"
        )
    print(f"Output               : {entry['output']}.png/.pdf")

    result = subprocess.run(
        make_command(entry, measure=False, vmin=vmin, vmax=vmax),
        env=env,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"find_centers_save.py failed while saving {entry['phase']}"
        )

print("============================================================")
print(f"Done. Files saved in: {output_dir}")
print("============================================================")
PY
