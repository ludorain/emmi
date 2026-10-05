#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

# ============================================================
# WHOLE-PICTURE EMMI ANNEALING PIPELINE
#
# Expected location:
#   emmi/whole_pictures_annealing_analysis.sh
#
# Expected input:
#   emmi/whole_pictures_annealing/A1/
#   emmi/whole_pictures_annealing/A2/
#   emmi/whole_pictures_annealing/B1/
#   emmi/whole_pictures_annealing/B2/
#
# Raw filenames:
#   A1_bef_ann_whole_run=..._data=denoised.tif
#   A1_bef_ann_whole_run=..._data=light.tif
#   A1_bef_ann_whole_run=..._data=denoised.png
#   A1_bef_ann_whole_run=..._data=light.png
#   A1_bef_ann_whole_run=..._overlay.png
#
#   A1_ann_75_5_whole_run=..._data=denoised.tif
#   ...
#
# Usage:
#   ./whole_pictures_annealing_analysis.sh
#   ./whole_pictures_annealing_analysis.sh A1
#
# With no argument all sensors are processed.
#
# IMPORTANT:
#   - NO luminosity macro is executed.
#   - Global hotspot IDs use only detected hotspot coordinates.
#   - All annealing phases are aligned to the before-annealing LIGHT image.
#   - The same measured geometrical transform is applied to LIGHT and DENOISED.
#   - Small whole-picture canvas mismatches are normalized before shift
#     measurement without image resizing/interpolation.
# ============================================================

# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$BASE_DIR/whole_pictures_annealing/pixels=36"

PROCESS_IMAGE="$BASE_DIR/manipulate_images/process_image.py"
TIF2TH2="$BASE_DIR/manipulate_images/tif2th2.py"

MEASURE_ROTATION="$BASE_DIR/tools_light_on/measure-rotation.py"
ROTATE_IMAGE="$BASE_DIR/tools_light_on/rotate-image.py"
MEASURE_SHIFT="$BASE_DIR/tools_light_on/measure-shift.py"
SHIFT_IMAGE="$BASE_DIR/tools_light_on/shift-image.py"

FIND_DEFECTS="$BASE_DIR/find_centers/find_defects_NOisolated_changeR.py"

PYTHON="${PYTHON:-python3}"

# Spatial tolerance for global-ID matching after alignment.
MATCH_RADIUS="${MATCH_RADIUS:-10.0}"

# Maximum allowed difference, in pixels per axis, between the rotated
# before-annealing reference canvas and a rotated moving-phase canvas.
# Small differences can arise because whole-picture acquisitions/stitching
# do not always contain exactly the same number of rows/columns.
# The pipeline normalizes such small mismatches WITHOUT resizing the image:
# larger canvases are center-cropped and smaller canvases are symmetrically
# padded. Larger mismatches remain fatal.
MAX_CANVAS_MISMATCH="${MAX_CANVAS_MISMATCH:-12}"

# If 1, every phase must contain exactly one data=light.tif and one
# data=denoised.tif. If 0, missing phases are left as "no_data".
STRICT_PHASES="${STRICT_PHASES:-1}"

# tif2th2.py was historically used with --error data=diffe.tif.
# Here no diffe image exists. We therefore try an input-only conversion of
# the final aligned denoised image. Failure is NON-FATAL because ROOT files
# are not needed for hotspot counting/global-ID assignment.
TRY_TIF2TH2="${TRY_TIF2TH2:0}"

PHASES=(
    "before_annealing"
    "annealing_75_5"
    "annealing_75_25"
    "annealing_100_5"
    "annealing_100_25"
    "annealing_125_5"
    "annealing_125_25"
    "annealing_150_5"
    "annealing_150_25"
)

ALL_SENSORS=(A1 A2 B1 B2)

# Sensor-specific state, reset before each sensor.
SENSOR=""
SENSOR_DIR=""
BASE_LIGHT_REFERENCE=""
ALIGNMENT_TABLE=""
MERGED_DIR=""

# ------------------------------------------------------------
# LOGGING / ERRORS
# ------------------------------------------------------------

log() {
    printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"
}

warn() {
    printf 'WARNING: %s\n' "$*" >&2
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

on_error() {
    local exit_code=$?
    local line_no=$1
    printf '\nERROR: pipeline stopped at line %s (exit code %s).\n' \
        "$line_no" "$exit_code" >&2
    exit "$exit_code"
}

trap 'on_error $LINENO' ERR

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------

SELECTED_SENSORS=()

if (( $# == 0 )); then
    SELECTED_SENSORS=("${ALL_SENSORS[@]}")
elif (( $# == 1 )); then
    case "$1" in
        A1|A2|B1|B2)
            SELECTED_SENSORS=("$1")
            ;;
        *)
            die "Invalid sensor '$1'. Allowed values: A1 A2 B1 B2."
            ;;
    esac
else
    cat >&2 <<EOF
Usage:
  $0
  $0 SENSOR

SENSOR:
  A1 | A2 | B1 | B2
EOF
    exit 1
fi

# ------------------------------------------------------------
# REQUIREMENTS
# ------------------------------------------------------------

check_requirements() {
    [[ -d "$DATA_DIR" ]] || die "Data directory not found: $DATA_DIR"

    command -v "$PYTHON" >/dev/null 2>&1 \
        || die "Python executable not found: $PYTHON"

    local required_files=(
        "$PROCESS_IMAGE"
        "$MEASURE_ROTATION"
        "$ROTATE_IMAGE"
        "$MEASURE_SHIFT"
        "$SHIFT_IMAGE"
        "$FIND_DEFECTS"
    )

    if [[ "$TRY_TIF2TH2" == "1" ]]; then
        required_files+=("$TIF2TH2")
    fi

    local f
    for f in "${required_files[@]}"; do
        [[ -f "$f" ]] || die "Required file not found: $f"
    done

    "$PYTHON" - <<'PY'
import numpy
import pandas
import tifffile
PY
}

# ------------------------------------------------------------
# PHASE / FILENAME HELPERS
# ------------------------------------------------------------

phase_prefix() {
    local sensor="$1"
    local phase="$2"

    case "$phase" in
        before_annealing)
            printf '%s\n' "${sensor}_bef_ann_whole"
            ;;
        annealing_75_5)
            printf '%s\n' "${sensor}_ann_75_5_whole"
            ;;
        annealing_75_25)
            printf '%s\n' "${sensor}_ann_75_25_whole"
            ;;
        annealing_100_5)
            printf '%s\n' "${sensor}_ann_100_5_whole"
            ;;
        annealing_100_25)
            printf '%s\n' "${sensor}_ann_100_25_whole"
            ;;
        annealing_125_5)
            printf '%s\n' "${sensor}_ann_125_5_whole"
            ;;
        annealing_125_25)
            printf '%s\n' "${sensor}_ann_125_25_whole"
            ;;
        annealing_150_5)
            printf '%s\n' "${sensor}_ann_150_5_whole"
            ;;
        annealing_150_25)
            printf '%s\n' "${sensor}_ann_150_25_whole"
            ;;
        *)
            die "Unknown phase: $phase"
            ;;
    esac
}

strip_tif_extension() {
    local name="$1"
    name="${name%.tif}"
    name="${name%.TIF}"
    name="${name%.tiff}"
    name="${name%.TIFF}"
    printf '%s\n' "$name"
}

extract_run_number_from_filename() {
    local path="$1"

    "$PYTHON" - "$path" <<'PY'
import os
import re
import sys

name = os.path.basename(sys.argv[1])
m = re.search(r"_run=([^_]+)", name)
if m is None:
    raise SystemExit(f"Cannot extract run number from filename: {name}")
print(m.group(1))
PY
}

get_unique_tif_by_kind() {
    local directory="$1"
    local kind="$2"

    local matches=()
    local f

    while IFS= read -r -d '' f; do
        matches+=("$f")
    done < <(
        find "$directory" \
            -maxdepth 1 \
            -type f \
            \( -iname '*.tif' -o -iname '*.tiff' \) \
            -name "*data=${kind}*" \
            -print0
    )

    (( ${#matches[@]} > 0 )) \
        || return 1

    if (( ${#matches[@]} != 1 )); then
        printf 'ERROR: expected exactly one data=%s TIF in %s, found %d:\n' \
            "$kind" "$directory" "${#matches[@]}" >&2
        printf '  %s\n' "${matches[@]}" >&2
        return 2
    fi

    printf '%s\n' "${matches[0]}"
}

phase_has_required_inputs() {
    local phase="$1"
    local originals_dir="$SENSOR_DIR/$phase/1originals"

    local light denoised

    if ! light="$(get_unique_tif_by_kind "$originals_dir" "light")"; then
        return 1
    fi

    if ! denoised="$(get_unique_tif_by_kind "$originals_dir" "denoised")"; then
        return 1
    fi

    [[ -n "$light" && -n "$denoised" ]]
}

# ------------------------------------------------------------
# STEP 1-3 - CREATE PHASE TREE AND MOVE RAW FILES
# ------------------------------------------------------------

organize_sensor_files() {
    log "Organizing raw files for sensor $SENSOR"

    [[ -d "$SENSOR_DIR" ]] \
        || die "Sensor directory not found: $SENSOR_DIR"

    local phase phase_dir originals_dir prefix
    local src target

    for phase in "${PHASES[@]}"; do
        phase_dir="$SENSOR_DIR/$phase"
        originals_dir="$phase_dir/1originals"

        mkdir -p \
            "$originals_dir" \
            "$phase_dir/2processed" \
            "$phase_dir/3rotated" \
            "$phase_dir/4coordinates" \
            "$phase_dir/5th2f"

        prefix="$(phase_prefix "$SENSOR" "$phase")"

        while IFS= read -r -d '' src; do
            target="$originals_dir/$(basename "$src")"

            if [[ -e "$target" ]]; then
                warn "Target already exists; leaving source untouched: $target"
                continue
            fi

            mv "$src" "$target"
            printf 'Moved: %s -> %s\n' "$(basename "$src")" "$phase/1originals/"
        done < <(
            find "$SENSOR_DIR" \
                -maxdepth 1 \
                -type f \
                -name "${prefix}_*" \
                -print0
        )
    done
}

# ------------------------------------------------------------
# OUTPUT-DIRECTORY RESET
# ------------------------------------------------------------

prepare_phase_output_dirs() {
    local phase="$1"
    local phase_dir="$SENSOR_DIR/$phase"

    rm -rf \
        "$phase_dir/2processed" \
        "$phase_dir/3rotated" \
        "$phase_dir/4coordinates" \
        "$phase_dir/5th2f" \
        "$phase_dir/.3rotated_pre_shift"

    mkdir -p \
        "$phase_dir/2processed" \
        "$phase_dir/3rotated" \
        "$phase_dir/4coordinates" \
        "$phase_dir/5th2f"
}

# ------------------------------------------------------------
# STEP 4 - IMAGE CLEANUP
# ------------------------------------------------------------

cleanup_images() {
    local phase="$1"
    local phase_dir="$SENSOR_DIR/$phase"
    local originals_dir="$phase_dir/1originals"
    local processed_dir="$phase_dir/2processed"

    log "Cleaning TIF images: $SENSOR / $phase"

    local n_files=0
    local input_file filename stem output_file

    while IFS= read -r -d '' input_file; do
        filename="$(basename "$input_file")"
        stem="$(strip_tif_extension "$filename")"
        output_file="$processed_dir/${stem}_processed.tif"

        "$PYTHON" "$PROCESS_IMAGE" \
            --input "$input_file" \
            --process remove_column_bias remove_hot_pixels remove_cold_pixels \
            --output "$output_file"

        ((n_files += 1))
    done < <(
        find "$originals_dir" \
            -maxdepth 1 \
            -type f \
            \( -iname '*.tif' -o -iname '*.tiff' \) \
            -print0
    )

    (( n_files > 0 )) \
        || die "No TIF files found in $originals_dir"
}

# ------------------------------------------------------------
# ROTATION HELPERS
# ------------------------------------------------------------

measure_rotation_angle() {
    local light_image="$1"
    local output angle

    output="$(
        "$PYTHON" "$MEASURE_ROTATION" --input "$light_image"
    )"

    printf '%s\n' "$output" >&2

    angle="$(
        printf '%s\n' "$output" |
        "$PYTHON" -c '
import re
import sys

text = sys.stdin.read()
m = re.search(
    r"rotation angle:\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
    text,
)
if m is None:
    raise SystemExit("Cannot parse rotation angle from measure-rotation.py output.")
print(m.group(1))
'
    )"

    printf '%s\n' "$angle"
}

measure_rotation_clip() {
    local input_image="$1"
    local rotated_image="$2"

    "$PYTHON" - "$input_image" "$rotated_image" <<'PY'
import sys
import tifffile

before = tifffile.imread(sys.argv[1])
after = tifffile.imread(sys.argv[2])

if before.ndim < 2 or after.ndim < 2:
    raise SystemExit("Rotation clipping check requires at least 2D images.")

ny0, nx0 = before.shape[-2:]
ny1, nx1 = after.shape[-2:]

dy = ny0 - ny1
dx = nx0 - nx1

if dy < 0 or dx < 0:
    clipy = 0.0
    clipx = 0.0
else:
    clipy = dy / 2.0
    clipx = dx / 2.0

print(f"{clipy:g} {clipx:g}")
PY
}

rotate_all_processed_images() {
    local processed_dir="$1"
    local destination_dir="$2"
    local angle="$3"

    mkdir -p "$destination_dir"

    local input_file filename stem output_file
    local n_rotated=0

    while IFS= read -r -d '' input_file; do
        filename="$(basename "$input_file")"
        stem="$(strip_tif_extension "$filename")"
        output_file="$destination_dir/${stem}_rotated.tif"

        "$PYTHON" "$ROTATE_IMAGE" \
            --input "$input_file" \
            --angle "$angle" \
            --output "$output_file"

        ((n_rotated += 1))
    done < <(
        find "$processed_dir" \
            -maxdepth 1 \
            -type f \
            -name '*_processed.tif' \
            -print0
    )

    (( n_rotated > 0 )) \
        || die "No processed images found in $processed_dir"
}

# ------------------------------------------------------------
# CANVAS-NORMALIZATION HELPER
# ------------------------------------------------------------
#
# Whole-picture acquisitions can differ by one or a few rows/columns even
# before rotation (for example 3046 versus 3047 columns). measure-shift.py
# requires identical array shapes, so every moving phase is normalized to the
# rotated before-annealing LIGHT canvas BEFORE measuring the translation.
#
# IMPORTANT:
#   - there is NO image rescaling/interpolation here;
#   - 1 input pixel remains 1 output pixel;
#   - if the moving image is larger, only border pixels are cropped;
#   - if it is smaller, border pixels are padded;
#   - LIGHT and DENOISED are normalized in exactly the same way because all
#     *_rotated.tif files in the phase pre-alignment directory are processed;
#   - a mismatch larger than MAX_CANVAS_MISMATCH is considered suspicious and
#     stops the pipeline.
#
# The crop/padding is centered. For an odd one-pixel excess, the extra pixel is
# removed from the bottom/right side. Any residual one-pixel origin difference
# is subsequently measured by measure-shift.py and applied to both images.
# ------------------------------------------------------------

normalize_rotated_canvas_to_reference() {
    local reference="$1"
    local directory="$2"

    "$PYTHON" - \
        "$reference" \
        "$directory" \
        "$MAX_CANVAS_MISMATCH" <<'PY_CANVAS'
from pathlib import Path
import sys

import numpy as np
import tifffile

reference_path = Path(sys.argv[1])
directory = Path(sys.argv[2])
max_mismatch = int(sys.argv[3])

reference = tifffile.imread(reference_path)

if reference.ndim != 2:
    raise SystemExit(
        f"Expected a 2D reference image, got shape {reference.shape}: "
        f"{reference_path}"
    )

target_y, target_x = reference.shape

files = sorted(directory.glob("*_rotated.tif"))
if not files:
    raise SystemExit(f"No rotated TIF images found in {directory}")

print("Canvas normalization:")
print(f"  reference = {reference_path.name}")
print(f"  target shape = {(target_y, target_x)}")
print(f"  maximum allowed mismatch = {max_mismatch} px/axis")

for path in files:
    image = tifffile.imread(path)

    if image.ndim != 2:
        raise SystemExit(
            f"Expected a 2D moving image, got shape {image.shape}: {path}"
        )

    original_shape = image.shape
    original_dtype = image.dtype
    ny, nx = original_shape

    # Positive values mean that the moving canvas is larger.
    dy = ny - target_y
    dx = nx - target_x

    if abs(dy) > max_mismatch or abs(dx) > max_mismatch:
        raise SystemExit(
            "Rotated canvas differs too much from the before-annealing "
            "reference:\n"
            f"  file:       {path}\n"
            f"  reference:  {(target_y, target_x)}\n"
            f"  moving:     {original_shape}\n"
            f"  difference: dy={dy}, dx={dx}\n"
            f"  allowed:    +/-{max_mismatch} px per axis\n"
            "This is probably not a harmless whole-picture/stitching size "
            "difference. Check the acquisition before increasing "
            "MAX_CANVAS_MISMATCH."
        )

    # --------------------------------------------------------
    # 1) CENTER-CROP dimensions that are larger than reference.
    # --------------------------------------------------------
    crop_top = crop_bottom = crop_left = crop_right = 0

    if image.shape[0] > target_y:
        excess = image.shape[0] - target_y
        crop_top = excess // 2
        crop_bottom = excess - crop_top
        image = image[crop_top:crop_top + target_y, :]

    if image.shape[1] > target_x:
        excess = image.shape[1] - target_x
        crop_left = excess // 2
        crop_right = excess - crop_left
        image = image[:, crop_left:crop_left + target_x]

    # --------------------------------------------------------
    # 2) SYMMETRICALLY PAD dimensions smaller than reference.
    # --------------------------------------------------------
    missing_y = target_y - image.shape[0]
    missing_x = target_x - image.shape[1]

    pad_top = max(missing_y, 0) // 2
    pad_bottom = max(missing_y, 0) - pad_top
    pad_left = max(missing_x, 0) // 2
    pad_right = max(missing_x, 0) - pad_left

    if missing_y > 0 or missing_x > 0:
        # Use a robust background-like value rather than zero so that a new
        # artificial black frame does not dominate the phase correlation.
        fill_value = np.median(image).item()
        image = np.pad(
            image,
            ((pad_top, pad_bottom), (pad_left, pad_right)),
            mode="constant",
            constant_values=fill_value,
        )

    if image.shape != (target_y, target_x):
        raise RuntimeError(
            f"Internal canvas-normalization error for {path}: "
            f"final shape={image.shape}, expected={(target_y, target_x)}"
        )

    # Preserve exactly the dtype of each moving TIFF.
    image = image.astype(original_dtype, copy=False)
    tifffile.imwrite(path, image)

    actions = []
    if any((crop_top, crop_bottom, crop_left, crop_right)):
        actions.append(
            "crop(top,bottom,left,right)="
            f"({crop_top},{crop_bottom},{crop_left},{crop_right})"
        )
    if any((pad_top, pad_bottom, pad_left, pad_right)):
        actions.append(
            "pad(top,bottom,left,right)="
            f"({pad_top},{pad_bottom},{pad_left},{pad_right})"
        )
    if not actions:
        actions.append("unchanged")

    print(
        f"  {path.name}: {original_shape} -> {image.shape}; "
        + "; ".join(actions)
    )

print("Canvas normalization completed.")
PY_CANVAS
}

# ------------------------------------------------------------
# SHIFT HELPERS
# ------------------------------------------------------------

assert_same_image_shape() {
    local reference="$1"
    local moving="$2"

    "$PYTHON" - "$reference" "$moving" <<'PY'
import sys
import tifffile

ref = tifffile.imread(sys.argv[1])
mov = tifffile.imread(sys.argv[2])

if ref.shape != mov.shape:
    raise SystemExit(
        "Reference and moving LIGHT images have different shapes after rotation:\n"
        f"  reference: {ref.shape}\n"
        f"  moving:    {mov.shape}\n"
        "measure-shift.py should only be run on geometrically compatible images."
    )
PY
}

measure_shift_values() {
    local reference="$1"
    local moving="$2"
    local output

    output="$(
        "$PYTHON" "$MEASURE_SHIFT" \
            --input "$reference" "$moving"
    )"

    printf '%s\n' "$output" >&2

    printf '%s\n' "$output" |
    "$PYTHON" -c '
import re
import sys

text = sys.stdin.read()
number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

m = re.search(
    rf"detected subpixel offset\s*\(y,\s*x\)\s*:\s*"
    rf"\[\s*({number})[\s,]+({number})\s*\]",
    text,
)

if m is None:
    raise SystemExit("Cannot parse shift from measure-shift.py output.")

print(m.group(1), m.group(2))
'
}

shift_all_rotated_images() {
    local source_dir="$1"
    local final_dir="$2"
    local shifty="$3"
    local shiftx="$4"

    mkdir -p "$final_dir"

    local input_file filename output_file
    local n_shifted=0

    while IFS= read -r -d '' input_file; do
        filename="$(basename "$input_file")"
        output_file="$final_dir/$filename"

        "$PYTHON" "$SHIFT_IMAGE" \
            --input "$input_file" \
            --shift "$shifty" "$shiftx" \
            --output "$output_file"

        ((n_shifted += 1))
    done < <(
        find "$source_dir" \
            -maxdepth 1 \
            -type f \
            -name '*_rotated.tif' \
            -print0
    )

    (( n_shifted > 0 )) \
        || die "No rotated images found in $source_dir"
}

# ------------------------------------------------------------
# STEP 5 - ALIGNMENT
# ------------------------------------------------------------

align_before_annealing() {
    local phase="before_annealing"
    local phase_dir="$SENSOR_DIR/$phase"
    local processed_dir="$phase_dir/2processed"
    local rotated_dir="$phase_dir/3rotated"

    log "Aligning reference phase: $SENSOR / $phase"

    local light_processed
    light_processed="$(get_unique_tif_by_kind "$processed_dir" "light")" \
        || die "Expected exactly one processed data=light image in $processed_dir"

    local angle
    angle="$(measure_rotation_angle "$light_processed")"

    rotate_all_processed_images \
        "$processed_dir" \
        "$rotated_dir" \
        "$angle"

    local light_name light_stem
    light_name="$(basename "$light_processed")"
    light_stem="$(strip_tif_extension "$light_name")"

    BASE_LIGHT_REFERENCE="$rotated_dir/${light_stem}_rotated.tif"

    [[ -f "$BASE_LIGHT_REFERENCE" ]] \
        || die "Reference LIGHT image was not created: $BASE_LIGHT_REFERENCE"

    local clipy clipx
    read -r clipy clipx < <(
        measure_rotation_clip "$light_processed" "$BASE_LIGHT_REFERENCE"
    )

    local raw_light
    raw_light="$(get_unique_tif_by_kind "$phase_dir/1originals" "light")"
    local run_number
    run_number="$(extract_run_number_from_filename "$raw_light")"

    printf '%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$SENSOR" "$phase" "$run_number" \
        "$angle" "$clipy" "$clipx" "0" "0" \
        >> "$ALIGNMENT_TABLE"

    printf 'Reference alignment:\n'
    printf '  sensor = %s\n' "$SENSOR"
    printf '  phase  = %s\n' "$phase"
    printf '  run    = %s\n' "$run_number"
    printf '  angle  = %s deg\n' "$angle"
    printf '  clipy  = %s px/side\n' "$clipy"
    printf '  clipx  = %s px/side\n' "$clipx"
    printf '  shifty = 0 px\n'
    printf '  shiftx = 0 px\n'
    printf '  LIGHT reference = %s\n' "$BASE_LIGHT_REFERENCE"
}

align_annealing_phase() {
    local phase="$1"
    local phase_dir="$SENSOR_DIR/$phase"
    local processed_dir="$phase_dir/2processed"
    local rotated_dir="$phase_dir/3rotated"
    local prealign_dir="$phase_dir/.3rotated_pre_shift"

    rm -rf "$prealign_dir"
    mkdir -p "$prealign_dir" "$rotated_dir"

    log "Aligning phase: $SENSOR / $phase"

    local light_processed
    light_processed="$(get_unique_tif_by_kind "$processed_dir" "light")" \
        || die "Expected exactly one processed data=light image in $processed_dir"

    local angle
    angle="$(measure_rotation_angle "$light_processed")"

    # First apply the phase-specific rotation to every processed TIF.
    rotate_all_processed_images \
        "$processed_dir" \
        "$prealign_dir" \
        "$angle"

    local light_name light_stem moving_light
    light_name="$(basename "$light_processed")"
    light_stem="$(strip_tif_extension "$light_name")"
    moving_light="$prealign_dir/${light_stem}_rotated.tif"

    [[ -f "$moving_light" ]] \
        || die "Rotated moving LIGHT image not found: $moving_light"

    local clipy clipx
    read -r clipy clipx < <(
        measure_rotation_clip "$light_processed" "$moving_light"
    )

    [[ -n "$BASE_LIGHT_REFERENCE" && -f "$BASE_LIGHT_REFERENCE" ]] \
        || die "Before-annealing LIGHT reference is unavailable."

    # Whole-picture runs may differ by one or a few pixels already at the raw
    # acquisition/stitching level. Normalize EVERY rotated image in this phase
    # (LIGHT and DENOISED) to the before-annealing reference canvas before
    # calling measure-shift.py. No resizing is performed.
    normalize_rotated_canvas_to_reference \
        "$BASE_LIGHT_REFERENCE" \
        "$prealign_dir"

    # The path is unchanged by normalization; verify that the LIGHT image now
    # has exactly the same shape as the reference.
    [[ -f "$moving_light" ]] \
        || die "Moving LIGHT disappeared during canvas normalization: $moving_light"

    assert_same_image_shape \
        "$BASE_LIGHT_REFERENCE" \
        "$moving_light"

    local shifty shiftx
    read -r shifty shiftx < <(
        measure_shift_values \
            "$BASE_LIGHT_REFERENCE" \
            "$moving_light"
    )

    # Apply the exact same translation to LIGHT and DENOISED.
    shift_all_rotated_images \
        "$prealign_dir" \
        "$rotated_dir" \
        "$shifty" \
        "$shiftx"

    rm -rf "$prealign_dir"

    local raw_light
    raw_light="$(get_unique_tif_by_kind "$phase_dir/1originals" "light")"
    local run_number
    run_number="$(extract_run_number_from_filename "$raw_light")"

    printf '%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$SENSOR" "$phase" "$run_number" \
        "$angle" "$clipy" "$clipx" "$shifty" "$shiftx" \
        >> "$ALIGNMENT_TABLE"

    printf 'Alignment parameters:\n'
    printf '  sensor = %s\n' "$SENSOR"
    printf '  phase  = %s\n' "$phase"
    printf '  run    = %s\n' "$run_number"
    printf '  angle  = %s deg\n' "$angle"
    printf '  clipy  = %s px/side\n' "$clipy"
    printf '  clipx  = %s px/side\n' "$clipx"
    printf '  shifty = %s px\n' "$shifty"
    printf '  shiftx = %s px\n' "$shiftx"
}

# ------------------------------------------------------------
# STEP 6 - HOTSPOT COORDINATES
# ------------------------------------------------------------

generate_phase_coordinates() {
    local phase="$1"
    local phase_dir="$SENSOR_DIR/$phase"
    local rotated_dir="$phase_dir/3rotated"
    local coordinates_dir="$phase_dir/4coordinates"

    local denoised_aligned
    denoised_aligned="$(get_unique_tif_by_kind "$rotated_dir" "denoised")" \
        || die "Expected exactly one aligned data=denoised image in $rotated_dir"

    local coordinate_file="$coordinates_dir/${SENSOR}_${phase}_coordinates.txt"

    log "Finding hotspots: $SENSOR / $phase"

    "$PYTHON" "$FIND_DEFECTS" \
        --input "$denoised_aligned" \
        --coordinates_root "$coordinate_file"

    [[ -s "$coordinate_file" ]] \
        || die "Coordinate file was not created or is empty: $coordinate_file"

    local n_hotspots
    n_hotspots="$(
        "$PYTHON" - "$coordinate_file" <<'PY'
import sys
from pathlib import Path

path = Path(sys.argv[1])
count = sum(1 for line in path.read_text().splitlines() if line.strip())
print(count)
PY
    )"

    printf 'Detected hotspots: %s / %s = %s\n' \
        "$SENSOR" "$phase" "$n_hotspots"
}

# ------------------------------------------------------------
# OPTIONAL: ALIGNED DENOISED TIF -> TH2F
# ------------------------------------------------------------

try_convert_denoised_to_th2f() {
    local phase="$1"

    [[ "$TRY_TIF2TH2" == "1" ]] || return 0

    local phase_dir="$SENSOR_DIR/$phase"
    local rotated_dir="$phase_dir/3rotated"
    local root_dir="$phase_dir/5th2f"

    local denoised_aligned
    denoised_aligned="$(get_unique_tif_by_kind "$rotated_dir" "denoised")" \
        || die "Expected exactly one aligned data=denoised image in $rotated_dir"

    local filename stem output_file
    filename="$(basename "$denoised_aligned")"
    stem="$(strip_tif_extension "$filename")"
    output_file="$root_dir/${stem}_th2f.root"

    log "Trying TIF -> TH2F conversion without an error map: $SENSOR / $phase"

    if "$PYTHON" "$TIF2TH2" \
        --input "$denoised_aligned" \
        --output "$output_file"
    then
        if [[ -s "$output_file" ]]; then
            printf 'Created optional TH2F: %s\n' "$output_file"
        else
            rm -f "$output_file"
            warn "tif2th2.py returned success but no non-empty ROOT file was produced."
        fi
    else
        rm -f "$output_file"
        warn "tif2th2.py could not convert without --error. Continuing: TH2F is not required by this pipeline."
    fi
}

# ------------------------------------------------------------
# PROCESS ONE PHASE
# ------------------------------------------------------------

process_phase() {
    local phase="$1"

    if ! phase_has_required_inputs "$phase"; then
        if [[ "$STRICT_PHASES" == "1" ]]; then
            die "Missing/ambiguous LIGHT or DENOISED TIF for $SENSOR / $phase"
        fi

        warn "Skipping $SENSOR / $phase because required input TIFs are missing."
        printf '%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "$SENSOR" "$phase" "//" "//" "//" "//" "//" "//" \
            >> "$ALIGNMENT_TABLE"
        return 0
    fi

    prepare_phase_output_dirs "$phase"
    cleanup_images "$phase"

    if [[ "$phase" == "before_annealing" ]]; then
        align_before_annealing
    else
        align_annealing_phase "$phase"
    fi

    generate_phase_coordinates "$phase"
    try_convert_denoised_to_th2f "$phase"
}

# ------------------------------------------------------------
# STEP 7 - GLOBAL HOTSPOT IDS + TEMPORAL DIAGNOSTICS
# ------------------------------------------------------------

merge_global_hotspot_ids() {
    mkdir -p "$MERGED_DIR"

    log "Global hotspot-ID assignment: $SENSOR"

    local phases_env
    phases_env="$(IFS='|'; printf '%s' "${PHASES[*]}")"

    SENSOR_ENV="$SENSOR" \
    SENSOR_DIR_ENV="$SENSOR_DIR" \
    MERGED_DIR_ENV="$MERGED_DIR" \
    MATCH_RADIUS_ENV="$MATCH_RADIUS" \
    PHASES_ENV="$phases_env" \
    "$PYTHON" <<'PY_GLOBAL'
from pathlib import Path
import math
import os
import re
import sys

import numpy as np
import pandas as pd

sensor = os.environ["SENSOR_ENV"]
sensor_dir = Path(os.environ["SENSOR_DIR_ENV"])
merged_dir = Path(os.environ["MERGED_DIR_ENV"])
match_radius = float(os.environ["MATCH_RADIUS_ENV"])
phases = os.environ["PHASES_ENV"].split("|")

merged_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def read_coordinates(path: Path) -> pd.DataFrame:
    rows = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue

            fields = [part.strip() for part in line.split(",")]
            if len(fields) < 2:
                raise ValueError(
                    f"Invalid coordinate line {line_number} in {path}: "
                    f"{raw.rstrip()}"
                )

            try:
                x = float(fields[0])
                y = float(fields[1])
                area = float(fields[2]) if len(fields) >= 3 else np.nan
                radius = float(fields[3]) if len(fields) >= 4 else np.nan
            except ValueError as exc:
                raise ValueError(
                    f"Non-numeric coordinate line {line_number} in {path}: "
                    f"{raw.rstrip()}"
                ) from exc

            rows.append({
                "local_spot": len(rows),
                "x": x,
                "y": y,
                "integration_area": area,
                "integration_radius": radius,
            })

    # Zero hotspots is a valid result in principle.
    return pd.DataFrame(
        rows,
        columns=[
            "local_spot",
            "x",
            "y",
            "integration_area",
            "integration_radius",
        ],
    )


def unique_coordinate_file(phase: str):
    coord_dir = sensor_dir / phase / "4coordinates"
    files = sorted(coord_dir.glob(f"{sensor}_{phase}_coordinates.txt"))

    if len(files) == 0:
        return None

    if len(files) != 1:
        raise RuntimeError(
            f"Expected exactly one coordinate file in {coord_dir}, "
            f"found {len(files)}."
        )

    return files[0]


def get_run_number(phase: str):
    originals = sensor_dir / phase / "1originals"
    light_files = sorted(originals.glob("*data=light.tif"))

    if not light_files:
        light_files = sorted(originals.glob("*data=light.tiff"))

    if len(light_files) != 1:
        return "//"

    m = re.search(r"_run=([^_]+)", light_files[0].name)
    return m.group(1) if m else "//"


def greedy_one_to_one(source: pd.DataFrame, target: pd.DataFrame, radius: float):
    """
    source columns: local_spot, x, y
    target columns: spot, x_ref, y_ref

    Returns:
        mapping[local_spot] = global_spot
        distance[local_spot] = distance
        candidate_count[local_spot] = number of candidates within radius
    """
    pairs = []
    candidate_count = {}

    for src in source.itertuples(index=False):
        n = 0
        for tgt in target.itertuples(index=False):
            d = math.hypot(
                float(src.x) - float(tgt.x_ref),
                float(src.y) - float(tgt.y_ref),
            )
            if d <= radius:
                pairs.append((d, int(src.local_spot), int(tgt.spot)))
                n += 1

        candidate_count[int(src.local_spot)] = n

    pairs.sort(key=lambda item: item[0])

    used_local = set()
    used_global = set()
    mapping = {}
    distances = {}

    for distance, local_spot, global_spot in pairs:
        if local_spot in used_local or global_spot in used_global:
            continue

        mapping[local_spot] = global_spot
        distances[local_spot] = distance
        used_local.add(local_spot)
        used_global.add(global_spot)

    return mapping, distances, candidate_count


# ------------------------------------------------------------
# Discover phase coordinate datasets
# ------------------------------------------------------------

phase_data = []

for phase_order, phase in enumerate(phases):
    coord_path = unique_coordinate_file(phase)

    if coord_path is None:
        phase_data.append({
            "phase": phase,
            "phase_order": phase_order,
            "available": False,
            "run_number": get_run_number(phase),
            "coord_path": None,
            "coords": None,
        })
        continue

    coords = read_coordinates(coord_path)

    phase_data.append({
        "phase": phase,
        "phase_order": phase_order,
        "available": True,
        "run_number": get_run_number(phase),
        "coord_path": coord_path,
        "coords": coords,
    })

if not phase_data[0]["available"]:
    raise RuntimeError(
        "before_annealing coordinates are required as the global-ID reference."
    )

# ------------------------------------------------------------
# Build global catalog
# ------------------------------------------------------------

catalog_rows = []
mapping_rows = []
phase_present = {}
phase_local_to_global = {}

reference = phase_data[0]
reference_coords = reference["coords"]

for row in reference_coords.itertuples(index=False):
    gid = len(catalog_rows)

    catalog_rows.append({
        "spot": gid,
        "x_ref": float(row.x),
        "y_ref": float(row.y),
        "first_seen_phase": reference["phase"],
        "first_seen_order": reference["phase_order"],
        "first_seen_run": reference["run_number"],
    })

    mapping_rows.append({
        "sensor": sensor,
        "phase": reference["phase"],
        "phase_order": reference["phase_order"],
        "run_number": reference["run_number"],
        "local_spot": int(row.local_spot),
        "spot": gid,
        "x": float(row.x),
        "y": float(row.y),
        "integration_area": row.integration_area,
        "integration_radius": row.integration_radius,
        "x_ref": float(row.x),
        "y_ref": float(row.y),
        "match_distance": 0.0,
        "match_type": "reference",
        "candidate_count": 1,
    })

phase_present[reference["phase"]] = set(range(len(catalog_rows)))
phase_local_to_global[reference["phase"]] = {
    int(row.local_spot): int(row.local_spot)
    for row in reference_coords.itertuples(index=False)
}

next_global_id = len(catalog_rows)

for d in phase_data[1:]:
    phase = d["phase"]

    if not d["available"]:
        phase_present[phase] = None
        phase_local_to_global[phase] = {}
        continue

    coords = d["coords"]

    catalog = pd.DataFrame(
        catalog_rows,
        columns=[
            "spot",
            "x_ref",
            "y_ref",
            "first_seen_phase",
            "first_seen_order",
            "first_seen_run",
        ],
    )

    mapping, distances, candidate_counts = greedy_one_to_one(
        source=coords,
        target=catalog,
        radius=match_radius,
    )

    match_types = {
        local_spot: "matched"
        for local_spot in mapping
    }

    # Every unmatched detected hotspot becomes a new global hotspot.
    for row in coords.itertuples(index=False):
        local_spot = int(row.local_spot)

        if local_spot in mapping:
            continue

        gid = next_global_id
        next_global_id += 1

        mapping[local_spot] = gid
        distances[local_spot] = np.nan
        candidate_counts.setdefault(local_spot, 0)
        match_types[local_spot] = "new"

        catalog_rows.append({
            "spot": gid,
            "x_ref": float(row.x),
            "y_ref": float(row.y),
            "first_seen_phase": phase,
            "first_seen_order": d["phase_order"],
            "first_seen_run": d["run_number"],
        })

        print(
            "NEW HOTSPOT: "
            f"sensor={sensor}, phase={phase}, "
            f"local_spot={local_spot}, global_spot={gid}, "
            f"x={float(row.x):.3f}, y={float(row.y):.3f}"
        )

    catalog_lookup = {
        int(row["spot"]): row
        for row in catalog_rows
    }

    for row in coords.itertuples(index=False):
        local_spot = int(row.local_spot)
        gid = int(mapping[local_spot])
        ref = catalog_lookup[gid]
        n_candidates = int(candidate_counts.get(local_spot, 0))

        mapping_rows.append({
            "sensor": sensor,
            "phase": phase,
            "phase_order": d["phase_order"],
            "run_number": d["run_number"],
            "local_spot": local_spot,
            "spot": gid,
            "x": float(row.x),
            "y": float(row.y),
            "integration_area": row.integration_area,
            "integration_radius": row.integration_radius,
            "x_ref": float(ref["x_ref"]),
            "y_ref": float(ref["y_ref"]),
            "match_distance": distances.get(local_spot, np.nan),
            "match_type": match_types[local_spot],
            "candidate_count": n_candidates,
        })

        if n_candidates > 1:
            print(
                "WARNING: multiple global-ID candidates within MATCH_RADIUS: "
                f"sensor={sensor}, phase={phase}, "
                f"local_spot={local_spot}, candidates={n_candidates}. "
                "Closest available match used.",
                file=sys.stderr,
            )

    phase_present[phase] = set(int(v) for v in mapping.values())
    phase_local_to_global[phase] = {
        int(k): int(v) for k, v in mapping.items()
    }

    # Per-phase detection-only global-ID table.
    per_phase_rows = []
    for row in coords.itertuples(index=False):
        local_spot = int(row.local_spot)
        gid = int(mapping[local_spot])

        per_phase_rows.append({
            "sensor": sensor,
            "phase": phase,
            "run_number": d["run_number"],
            "local_spot": local_spot,
            "spot": gid,
            "x": float(row.x),
            "y": float(row.y),
            "integration_area": row.integration_area,
            "integration_radius": row.integration_radius,
            "match_type": match_types[local_spot],
            "match_distance": distances.get(local_spot, np.nan),
            "candidate_count": int(candidate_counts.get(local_spot, 0)),
        })

    per_phase_output = (
        merged_dir / f"{sensor}_{phase}_global_ID.csv"
    )
    pd.DataFrame(per_phase_rows).sort_values(
        "spot", kind="stable"
    ).to_csv(per_phase_output, index=False)
    print(f"Created: {per_phase_output}")

# Also create the before-annealing per-phase file.
before_rows = [
    row for row in mapping_rows
    if row["phase"] == "before_annealing"
]
before_output = merged_dir / f"{sensor}_before_annealing_global_ID.csv"
pd.DataFrame(before_rows).drop(
    columns=["phase_order"], errors="ignore"
).sort_values("spot", kind="stable").to_csv(before_output, index=False)
print(f"Created: {before_output}")

all_global_spots = sorted(int(row["spot"]) for row in catalog_rows)
catalog = pd.DataFrame(catalog_rows).sort_values("spot").reset_index(drop=True)
mapping_df = pd.DataFrame(mapping_rows)

# ------------------------------------------------------------
# Presence/status table and reappearance source
# ------------------------------------------------------------

first_seen_order = {
    int(row.spot): int(row.first_seen_order)
    for row in catalog.itertuples(index=False)
}

presence_rows = []

for gid in all_global_spots:
    last_detected_phase = None
    last_detected_order = None

    for i, d in enumerate(phase_data):
        phase = d["phase"]
        available = d["available"]
        present_set = phase_present.get(phase)

        if not available or present_set is None:
            presence_rows.append({
                "sensor": sensor,
                "spot": gid,
                "phase": phase,
                "phase_order": i,
                "run_number": d["run_number"],
                "detected": "//",
                "status": "no_data",
                "reappeared_from_phase": "//",
                "local_spot": "//",
            })
            continue

        detected = gid in present_set
        inverse = {
            global_spot: local_spot
            for local_spot, global_spot
            in phase_local_to_global[phase].items()
        }
        local_spot = inverse.get(gid, np.nan)

        if i == 0:
            if detected:
                status = "reference"
            else:
                status = "not_yet_detected"
            reappeared_from = "//"

        elif not phase_data[i - 1]["available"]:
            # There is no measured immediately preceding phase, therefore
            # disappearance/reappearance cannot be established.
            if detected:
                status = "detected_after_no_data"
                if last_detected_phase is None:
                    status = "first_detected_after_no_data"
            else:
                status = "absent_after_no_data"
            reappeared_from = "//"

        else:
            prev_phase = phase_data[i - 1]["phase"]
            prev_present_set = phase_present.get(prev_phase)
            was_detected_previous = (
                prev_present_set is not None and gid in prev_present_set
            )

            if i < first_seen_order[gid]:
                status = "not_yet_detected"
                reappeared_from = "//"

            elif i == first_seen_order[gid]:
                status = "appeared"
                reappeared_from = "//"

            elif detected and was_detected_previous:
                status = "present"
                reappeared_from = "//"

            elif detected and not was_detected_previous:
                # The hotspot existed previously, disappeared for >=1 phase,
                # and is detected again now.
                status = "reappeared"
                reappeared_from = (
                    last_detected_phase
                    if last_detected_phase is not None
                    else "//"
                )

            elif (not detected) and was_detected_previous:
                status = "disappeared"
                reappeared_from = "//"

            else:
                status = "absent"
                reappeared_from = "//"

        presence_rows.append({
            "sensor": sensor,
            "spot": gid,
            "phase": phase,
            "phase_order": i,
            "run_number": d["run_number"],
            "detected": bool(detected),
            "status": status,
            "reappeared_from_phase": reappeared_from,
            "local_spot": local_spot,
        })

        if detected:
            last_detected_phase = phase
            last_detected_order = i

presence = pd.DataFrame(presence_rows)

# ------------------------------------------------------------
# Last-seen information for catalog
# ------------------------------------------------------------

last_seen_phase = {}
last_seen_run = {}

for gid in all_global_spots:
    seen = presence[
        (presence["spot"] == gid) &
        (presence["detected"] == True)  # noqa: E712
    ]

    if seen.empty:
        last_seen_phase[gid] = "//"
        last_seen_run[gid] = "//"
    else:
        last = seen.iloc[-1]
        last_seen_phase[gid] = last["phase"]
        last_seen_run[gid] = last["run_number"]

catalog["last_seen_phase"] = catalog["spot"].map(last_seen_phase)
catalog["last_seen_run"] = catalog["spot"].map(last_seen_run)

# ------------------------------------------------------------
# Wide, human-readable phase summary requested by the user
# ------------------------------------------------------------

def reappearance_column_name(phase: str) -> str:
    return "reappeared_from_" + phase

summary_rows = []

for i, d in enumerate(phase_data):
    phase = d["phase"]

    row = {
        "sensor": sensor,
        "phase": phase,
        "run_number": d["run_number"],
    }

    if not d["available"]:
        row["counted_hotspots"] = "//"
        row["disappeared"] = "//"
        row["reappeared_total"] = "//"

        for source_phase in phases:
            row[reappearance_column_name(source_phase)] = "//"

        summary_rows.append(row)
        continue

    pp = presence[presence["phase"] == phase]

    row["counted_hotspots"] = int(
        sum(value is True or value == True for value in pp["detected"])  # noqa: E712
    )

    # There is no "previous phase" for before_annealing.
    # If the previous phase has no data, transition counts are undefined.
    transitions_defined = (
        i > 0 and phase_data[i - 1]["available"]
    )

    if not transitions_defined:
        row["disappeared"] = "//"
        row["reappeared_total"] = "//"
    else:
        row["disappeared"] = int((pp["status"] == "disappeared").sum())
        row["reappeared_total"] = int((pp["status"] == "reappeared").sum())

    for source_index, source_phase in enumerate(phases):
        col = reappearance_column_name(source_phase)

        # A hotspot cannot "reappear from" its current or a future phase.
        # The user requested these meaningless cells to contain //.
        if source_index >= i:
            row[col] = "//"
            continue

        if not transitions_defined:
            row[col] = "//"
            continue

        row[col] = int(
            (
                (pp["status"] == "reappeared") &
                (pp["reappeared_from_phase"] == source_phase)
            ).sum()
        )

    summary_rows.append(row)

summary = pd.DataFrame(summary_rows)

# ------------------------------------------------------------
# Save products
# ------------------------------------------------------------

mapping_output = merged_dir / f"{sensor}_global_mapping.csv"
mapping_df.drop(columns=["phase_order"], errors="ignore").sort_values(
    ["spot", "phase"], kind="stable"
).to_csv(mapping_output, index=False)
print(f"Created: {mapping_output}")

presence_output = merged_dir / f"{sensor}_hotspot_presence.csv"
presence.drop(columns=["phase_order"], errors="ignore").to_csv(
    presence_output, index=False
)
print(f"Created: {presence_output}")

catalog_output = merged_dir / f"{sensor}_global_catalog.csv"
catalog.drop(columns=["first_seen_order"], errors="ignore").to_csv(
    catalog_output, index=False
)
print(f"Created: {catalog_output}")

summary_output = merged_dir / f"{sensor}_hotspot_phase_summary.csv"
summary.to_csv(summary_output, index=False)
print(f"Created: {summary_output}")

print("")
print("============================================")
print(f"GLOBAL HOTSPOT REPORT - {sensor}")
print(f"Global hotspots: {len(all_global_spots)}")
print("")

for row in summary.itertuples(index=False):
    print(
        f"{row.phase}: "
        f"counted={row.counted_hotspots}, "
        f"disappeared={row.disappeared}, "
        f"reappeared={row.reappeared_total}"
    )

print("")
print(f"Summary:  {summary_output}")
print(f"Presence: {presence_output}")
print(f"Mapping:  {mapping_output}")
print(f"Catalog:  {catalog_output}")
print("============================================")
PY_GLOBAL
}

# ------------------------------------------------------------
# PROCESS ONE SENSOR
# ------------------------------------------------------------

process_sensor() {
    SENSOR="$1"
    SENSOR_DIR="$DATA_DIR/$SENSOR"
    MERGED_DIR="$SENSOR_DIR/merged_files"
    ALIGNMENT_TABLE="$MERGED_DIR/${SENSOR}_alignment.csv"
    BASE_LIGHT_REFERENCE=""

    printf '\n'
    printf '============================================================\n'
    printf 'PROCESSING SENSOR %s\n' "$SENSOR"
    printf 'Sensor directory: %s\n' "$SENSOR_DIR"
    printf 'Match radius:     %s px\n' "$MATCH_RADIUS"
    printf 'Canvas tolerance: %s px/axis\n' "$MAX_CANVAS_MISMATCH"
    printf '============================================================\n'

    [[ -d "$SENSOR_DIR" ]] \
        || die "Sensor directory not found: $SENSOR_DIR"

    organize_sensor_files

    mkdir -p "$MERGED_DIR"
    printf '%s\n' \
        "sensor,phase,run_number,angle,clipy,clipx,shifty,shiftx" \
        > "$ALIGNMENT_TABLE"

    local phase
    for phase in "${PHASES[@]}"; do
        process_phase "$phase"
    done

    merge_global_hotspot_ids

    printf '\nCompleted sensor %s\n' "$SENSOR"
    printf 'Alignment table: %s\n' "$ALIGNMENT_TABLE"
    printf 'Merged outputs:  %s\n' "$MERGED_DIR"
}

# ------------------------------------------------------------
# COMBINE SENSOR SUMMARIES
# ------------------------------------------------------------

build_all_sensor_summary() {
    local output="$DATA_DIR/hotspot_phase_summary_all_sensors.csv"

    DATA_DIR_ENV="$DATA_DIR" \
    OUTPUT_ENV="$output" \
    "$PYTHON" <<'PY'
from pathlib import Path
import os
import pandas as pd

data_dir = Path(os.environ["DATA_DIR_ENV"])
output = Path(os.environ["OUTPUT_ENV"])

frames = []

for sensor in ("A1", "A2", "B1", "B2"):
    path = (
        data_dir / sensor / "merged_files" /
        f"{sensor}_hotspot_phase_summary.csv"
    )
    if path.is_file():
        frames.append(pd.read_csv(path, dtype=str, keep_default_na=False))

if not frames:
    raise SystemExit("No per-sensor hotspot phase summaries were found.")

combined = pd.concat(frames, ignore_index=True, sort=False)
combined.to_csv(output, index=False)
print(f"Created combined summary: {output}")
PY
}

# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

main() {
    check_requirements

    printf '\n'
    printf '============================================================\n'
    printf 'WHOLE-PICTURE EMMI ANNEALING PIPELINE\n'
    printf 'BASE_DIR: %s\n' "$BASE_DIR"
    printf 'DATA_DIR: %s\n' "$DATA_DIR"
    printf 'Sensors:  %s\n' "${SELECTED_SENSORS[*]}"
    printf '============================================================\n'

    local sensor
    for sensor in "${SELECTED_SENSORS[@]}"; do
        process_sensor "$sensor"
    done

    build_all_sensor_summary

    printf '\n'
    printf '============================================================\n'
    printf 'PIPELINE COMPLETED SUCCESSFULLY\n'
    printf 'Combined summary:\n'
    printf '  %s\n' "$DATA_DIR/hotspot_phase_summary_all_sensors_pixel=36.csv"
    printf '============================================================\n'
}

main "$@"
