#!/usr/bin/env bash

# ============================================================
# DISPLAY ROTATION / SHIFT VALIDATION FOR ONE GLOBAL HOTSPOT
#
# Usage:
#   ./display_rotation_shift_pipeline.sh SENSOR CONSTANT GLOBAL_SPOT_ID PHASE
#
# Examples:
#   ./display_rotation_shift_pipeline.sh A1 T 31 annealing_T=75_h=5
#   ./display_rotation_shift_pipeline.sh B1 v 12 annealing_T=125_h=25
#
# The script must be located inside the main emmi directory.
# display-rotation-shift.py must be located in:
#   emmi/tools_light_on/display-rotation-shift.py
#
# Input data are read only from DATA_irradiated_isolated_changeR.
# No image-processing step is repeated.
# ============================================================

set -Eeuo pipefail
shopt -s nullglob

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${DATA_DIR:-$BASE_DIR/DATA_irradiated_isolated_changeR}"
DISPLAY_PY="${DISPLAY_PY:-$BASE_DIR/tools_light_on/display-rotation-shift.py}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$BASE_DIR/display_rotation_shift_results}"
PYTHON="${PYTHON:-python3}"

log() {
    printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

on_error() {
    local exit_code=$?
    local line_no=$1
    printf '\nERROR: display pipeline stopped at line %s (exit code %s).\n' \
        "$line_no" "$exit_code" >&2
    exit "$exit_code"
}
trap 'on_error $LINENO' ERR

if [[ $# -ne 4 ]]; then
    cat >&2 <<EOF
Usage:
  $0 SENSOR CONSTANT GLOBAL_SPOT_ID PHASE

Examples:
  $0 A1 T 31 annealing_T=75_h=5
  $0 B1 v 12 annealing_T=125_h=25
EOF
    exit 1
fi

SENSOR="$1"
CONSTANT="$2"
GLOBAL_SPOT="$3"
PHASE="$4"

case "$SENSOR" in
    A1|A2|B1|B2) ;;
    *) die "Invalid sensor '$SENSOR'. Allowed: A1 A2 B1 B2." ;;
esac

case "$CONSTANT" in
    T|v) ;;
    *) die "Invalid constant '$CONSTANT'. Allowed: T or v." ;;
esac

[[ "$GLOBAL_SPOT" =~ ^[0-9]+$ ]] \
    || die "GLOBAL_SPOT_ID must be a non-negative integer."

command -v "$PYTHON" >/dev/null 2>&1 \
    || die "Python executable not found: $PYTHON"
[[ -d "$DATA_DIR" ]] || die "Data directory not found: $DATA_DIR"
[[ -f "$DISPLAY_PY" ]] || die "Display program not found: $DISPLAY_PY"

find_unique_run() {
    local phase_dir="$1"
    local matches=()
    local d

    while IFS= read -r -d '' d; do
        matches+=("$d")
    done < <(
        find "$phase_dir" -mindepth 1 -maxdepth 1 -type d \
            -name "${SENSOR}_${CONSTANT}=*_run=*" -print0
    )

    (( ${#matches[@]} > 0 )) \
        || die "No ${SENSOR}_${CONSTANT}=*_run=* directory found in $phase_dir"
    (( ${#matches[@]} == 1 )) \
        || die "More than one matching run found in $phase_dir"

    printf '%s\n' "${matches[0]}"
}

REFERENCE_PHASE_DIR="$DATA_DIR/before_annealing"
TARGET_PHASE_DIR="$DATA_DIR/$PHASE"
[[ -d "$REFERENCE_PHASE_DIR" ]] || die "Missing before_annealing directory."
[[ -d "$TARGET_PHASE_DIR" ]] || die "Phase directory not found: $TARGET_PHASE_DIR"

REFERENCE_RUN="$(find_unique_run "$REFERENCE_PHASE_DIR")"
TARGET_RUN="$(find_unique_run "$TARGET_PHASE_DIR")"
REFERENCE_RUN_NAME="$(basename "$REFERENCE_RUN")"
TARGET_RUN_NAME="$(basename "$TARGET_RUN")"

if [[ ! "$REFERENCE_RUN_NAME" =~ ^(${SENSOR}_${CONSTANT}=([^_]+))_run=(.+)$ ]]; then
    die "Cannot parse reference run name: $REFERENCE_RUN_NAME"
fi
RUN_PREFIX="${BASH_REMATCH[1]}"
FIXED_VALUE="${BASH_REMATCH[2]}"
REFERENCE_RUN_NUMBER="${BASH_REMATCH[3]}"

if [[ ! "$TARGET_RUN_NAME" =~ ^(${SENSOR}_${CONSTANT}=([^_]+))_run=(.+)$ ]]; then
    die "Cannot parse target run name: $TARGET_RUN_NAME"
fi
TARGET_FIXED_VALUE="${BASH_REMATCH[2]}"
TARGET_RUN_NUMBER="${BASH_REMATCH[3]}"

[[ "$TARGET_FIXED_VALUE" == "$FIXED_VALUE" ]] || die \
    "Fixed ${CONSTANT} differs between before_annealing (${FIXED_VALUE}) and ${PHASE} (${TARGET_FIXED_VALUE})."

MERGED_DIR="$DATA_DIR/merged_files/$RUN_PREFIX"
GEOMETRY_PLAN="$MERGED_DIR/${SENSOR}_${CONSTANT}_geometry_plan.csv"
[[ -f "$GEOMETRY_PLAN" ]] || die "Geometry plan not found: $GEOMETRY_PLAN"

# ------------------------------------------------------------
# GLOBAL HOTSPOT REFERENCE COORDINATES
# ------------------------------------------------------------

GEOMETRY_INFO="$($PYTHON - "$GEOMETRY_PLAN" "$GLOBAL_SPOT" "$PHASE" <<'PY'
import csv
import sys

path, spot_text, target_phase = sys.argv[1:]
spot = int(spot_text)

rows = []
with open(path, newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        try:
            row_spot = int(float(row["spot"]))
        except (KeyError, ValueError, TypeError):
            continue
        if row_spot == spot:
            rows.append(row)

if not rows:
    raise SystemExit(f"Global hotspot {spot} is not present in {path}")

ref = next((r for r in rows if r.get("phase") == "before_annealing"), None)
target = next((r for r in rows if r.get("phase") == target_phase), None)

if ref is None:
    raise SystemExit(f"No before_annealing geometry row for global hotspot {spot}")
if target is None:
    raise SystemExit(
        f"No geometry row for phase {target_phase!r}, global hotspot {spot}"
    )

print(
    ref["x"],
    ref["y"],
    ref.get("detection", "unknown"),
    target.get("detection", "unknown"),
    sep="\t",
)
PY
)"

IFS=$'\t' read -r X_REF Y_REF REF_DETECTION TARGET_DETECTION <<< "$GEOMETRY_INFO"

[[ -n "$X_REF" && -n "$Y_REF" ]] \
    || die "Could not read reference coordinates for global hotspot $GLOBAL_SPOT."

# ------------------------------------------------------------
# FIND MATCHING data=diff IMAGES
# ------------------------------------------------------------
# Same operating-point rule as source finding:
#   T fixed -> second-highest v
#   v fixed -> second-highest T

IMAGE_INFO="$($PYTHON - \
    "$REFERENCE_RUN/3rotated" \
    "$TARGET_RUN/2processed" \
    "$TARGET_RUN/3rotated" \
    "$CONSTANT" <<'PY'
from pathlib import Path
import re
import sys

ref_dir = Path(sys.argv[1])
before_dir = Path(sys.argv[2])
after_dir = Path(sys.argv[3])
constant = sys.argv[4]

varying = "v" if constant == "T" else "T"
number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
pattern = re.compile(rf"(?:^|_){re.escape(varying)}=({number})")


def diff_files(directory: Path, suffix_hint: str):
    if not directory.is_dir():
        raise SystemExit(f"Missing image directory: {directory}")

    entries = []
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        if path.suffix.lower() not in {".tif", ".tiff"}:
            continue
        if "data=diff" not in path.name or "data=diffe" in path.name:
            continue
        if suffix_hint and suffix_hint not in path.name:
            continue

        match = pattern.search(path.name)
        if match is not None:
            entries.append((float(match.group(1)), path))

    return entries


ref_entries = diff_files(ref_dir, "_processed_rotated.tif")
if not ref_entries:
    raise SystemExit(f"No reference data=diff images found in {ref_dir}")

values = sorted({value for value, _path in ref_entries}, reverse=True)
if len(values) < 2:
    raise SystemExit(
        f"Cannot select second-highest {varying}: only {len(values)} "
        f"distinct value(s) in {ref_dir}"
    )

selected = values[1]


def one_at_value(entries, directory: Path):
    matches = [
        path
        for value, path in entries
        if abs(value - selected) < 1e-9
    ]

    if len(matches) != 1:
        raise SystemExit(
            f"Expected exactly one data=diff image at {varying}={selected:g} "
            f"in {directory}; found {len(matches)}"
        )

    return matches[0]


reference_image = one_at_value(ref_entries, ref_dir)
target_before_image = one_at_value(
    diff_files(before_dir, "_processed.tif"),
    before_dir,
)
target_after_image = one_at_value(
    diff_files(after_dir, "_processed_rotated.tif"),
    after_dir,
)

print(
    reference_image,
    target_before_image,
    target_after_image,
    f"{varying}={selected:g}",
    sep="\t",
)
PY
)"

IFS=$'\t' read -r REFERENCE_IMAGE TARGET_BEFORE_IMAGE TARGET_AFTER_IMAGE SCAN_LABEL <<< "$IMAGE_INFO"

[[ -f "$REFERENCE_IMAGE" ]] || die "Reference image not found: $REFERENCE_IMAGE"
[[ -f "$TARGET_BEFORE_IMAGE" ]] || die "Target pre-alignment image not found: $TARGET_BEFORE_IMAGE"
[[ -f "$TARGET_AFTER_IMAGE" ]] || die "Target aligned image not found: $TARGET_AFTER_IMAGE"

# ------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------

OUTPUT_DIR="$OUTPUT_ROOT/$PHASE/$TARGET_RUN_NAME/spot${GLOBAL_SPOT}"
mkdir -p "$OUTPUT_DIR"
find "$OUTPUT_DIR" -maxdepth 1 -type f -name '*.png' -delete

printf '\n============================================================\n'
printf 'ROTATION / SHIFT DISPLAY\n'
printf 'Sensor:             %s\n' "$SENSOR"
printf 'Constant:           %s=%s\n' "$CONSTANT" "$FIXED_VALUE"
printf 'Global hotspot:     %s\n' "$GLOBAL_SPOT"
printf 'Phase:              %s\n' "$PHASE"
printf 'Reference run:      %s\n' "$REFERENCE_RUN_NUMBER"
printf 'Target run:         %s\n' "$TARGET_RUN_NUMBER"
printf 'Reference center:   x=%s, y=%s\n' "$X_REF" "$Y_REF"
printf 'Reference detected: %s\n' "$REF_DETECTION"
printf 'Target detected:    %s\n' "$TARGET_DETECTION"
printf 'Displayed scan:     %s\n' "$SCAN_LABEL"
printf 'Display program:    %s\n' "$DISPLAY_PY"
printf 'Output:             %s\n' "$OUTPUT_DIR"
printf '============================================================\n'

if [[ "$REF_DETECTION" != "present" ]]; then
    printf 'WARNING: hotspot %s is not detected in before_annealing; the reference crosshair uses propagated geometry.\n' \
        "$GLOBAL_SPOT" >&2
fi

if [[ "$TARGET_DETECTION" != "present" ]]; then
    printf 'WARNING: hotspot %s is not detected in phase %s; a visible target hotspot may be absent.\n' \
        "$GLOBAL_SPOT" "$PHASE" >&2
fi

log "Creating 4-panel 50x50 px comparison display"

"$PYTHON" "$DISPLAY_PY" \
    --reference "$REFERENCE_IMAGE" \
    --target-before "$TARGET_BEFORE_IMAGE" \
    --target-after "$TARGET_AFTER_IMAGE" \
    --x "$X_REF" \
    --y "$Y_REF" \
    --sensor "$SENSOR" \
    --constant-label "${CONSTANT}=${FIXED_VALUE}" \
    --phase "$PHASE" \
    --spot "$GLOBAL_SPOT" \
    --scan-label "$SCAN_LABEL" \
    --output-dir "$OUTPUT_DIR" \
    --crop-size 60

printf '\nCreated comparison figure in:\n  %s\n' "$OUTPUT_DIR"
