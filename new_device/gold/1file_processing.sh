#!/usr/bin/env bash

# =============================================================================
# 1file_processing.sh
#
# Reduced-sample EMMI pipeline for:
#   emmi/new_device/gold/DATA/A1_T=20_run=20250826-221930
#   emmi/new_device/gold/DATA/A1_v=7_run=20250826-035836
#   emmi/new_device/gold/DATA/A1_v=9_run=20250825-093810
#
# Workflow
# --------
#  1. prepare each dataset locally: clean data=diff and data=light only;
#     data=diffe is NEVER cleaned, rotated or shifted.
#  2. find hotspots in each dataset native frame using the cleaned data=diff
#     image at the maximum varying operating parameter.
#  3. build TH2F locally from cleaned data=diff + ORIGINAL data=diffe.
#  4. compute luminosity/statistical error at R=20,16,24 using the same local
#     hotspot centre for every operating point of that dataset.
#  5. compute deltaL=max(|L24-L20|,|L20-L16|) while IDs are still local.
#  6. only after all measurements are complete, align the light images and
#     transform hotspot coordinates into a common frame for spatial matching.
#  7. assign global IDs. v9 seeds the numbering. Missing hotspots are never
#     force-measured.
#
# Measurement frame and matching frame are deliberately separate.
# =============================================================================

set -Eeuo pipefail
shopt -s nullglob

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EMMI_DIR="${EMMI_DIR:-$(cd "$BASE_DIR/../.." && pwd)}"
DATA_DIR="${DATA_DIR:-$BASE_DIR/DATA}"

SENSOR="A1"

RUN_T20="${RUN_T20:-$DATA_DIR/A1_T=20_run=20250826-221930}"
RUN_V7="${RUN_V7:-$DATA_DIR/A1_v=7_run=20250826-035836}"
RUN_V9="${RUN_V9:-$DATA_DIR/A1_v=9_run=20250825-093810}"

# Primary geometrical reference for inter-run alignment and final global-ID seed.
# This reference is NOT used as the luminosity-integration centre for other runs.
PRIMARY_RUN="$RUN_V9"
PRIMARY_KEY="v9"

PROCESS_IMAGE="$EMMI_DIR/manipulate_images/process_image.py"
TIF2TH2="$EMMI_DIR/manipulate_images/tif2th2.py"

MEASURE_ROTATION="$EMMI_DIR/tools_light_on/measure-rotation.py"
ROTATE_IMAGE="$EMMI_DIR/tools_light_on/rotate-image.py"
MEASURE_SHIFT="$EMMI_DIR/tools_light_on/measure-shift.py"
SHIFT_IMAGE="$EMMI_DIR/tools_light_on/shift-image.py"

# Same source finder used by the full automatic_file_processing pipeline.
# Its detected centres are authoritative for source finding; luminosity apertures
# are overridden below to the fixed nominal radius R=20.
FIND_DEFECTS="${FIND_DEFECTS:-$EMMI_DIR/find_centers/find_defects_isolated_changeR.py}"

SPOT_LUM_DIR="$EMMI_DIR/spot_luminosity"
SPOT_LUM_MACRO="${SPOT_LUM_MACRO:-spot_luminosity_sum_irradiated.C}"

PYTHON="${PYTHON:-python3}"

# Spatial tolerance used ONLY in the final global-ID assignment, after all
# luminosity and systematic-uncertainty calculations are complete.
MATCH_RADIUS="${MATCH_RADIUS:-10.0}"

# Tolerance only for reconnecting ROOT output rows to coordinates supplied to ROOT.
COORD_MATCH_RADIUS="${COORD_MATCH_RADIUS:-1.0}"

# Nominal and systematic-variation radii [pixel].
R_NOMINAL="20"
R_LOW="16"
R_HIGH="24"

# A1 breakdown voltage. Used only for the T=20 voltage scan to add v_fin.
VBD_A1="${VBD_A1:-51.3}"

# Radius-specific output trees. Nominal R=20 stays in DATA.
R16_DATA_DIR="${R16_DATA_DIR:-$BASE_DIR/DATA_R=16}"
R24_DATA_DIR="${R24_DATA_DIR:-$BASE_DIR/DATA_R=24}"

MERGED_DIR="$DATA_DIR/merged_files"
MANIFEST="$MERGED_DIR/A1_dataset_manifest.csv"
ALIGNMENT_TABLE="$MERGED_DIR/A1_alignment.csv"

# Filled after the primary run is rotated.
BASE_LIGHT_REFERENCE=""
LAST_SOURCE_COORDS=""

# -----------------------------------------------------------------------------
# LOGGING / ERRORS
# -----------------------------------------------------------------------------

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
    local code=$?
    local line_no=$1
    printf '\nERROR: 1file_processing.sh stopped at line %s (exit code %s).\n' \
        "$line_no" "$code" >&2
    exit "$code"
}

trap 'on_error $LINENO' ERR

# -----------------------------------------------------------------------------
# REQUIREMENTS
# -----------------------------------------------------------------------------

check_requirements() {
    command -v "$PYTHON" >/dev/null 2>&1 \
        || die "Python executable not found: $PYTHON"
    command -v root >/dev/null 2>&1 \
        || die "ROOT executable 'root' not found in PATH."

    local required_files=(
        "$PROCESS_IMAGE"
        "$TIF2TH2"
        "$MEASURE_ROTATION"
        "$ROTATE_IMAGE"
        "$MEASURE_SHIFT"
        "$SHIFT_IMAGE"
        "$FIND_DEFECTS"
        "$SPOT_LUM_DIR/$SPOT_LUM_MACRO"
    )

    local f
    for f in "${required_files[@]}"; do
        [[ -f "$f" ]] || die "Required file not found: $f"
    done

    local run
    for run in "$RUN_T20" "$RUN_V7" "$RUN_V9"; do
        [[ -d "$run/1originals" ]] || die "Missing input directory: $run/1originals"
    done
}

# -----------------------------------------------------------------------------
# GENERIC HELPERS
# -----------------------------------------------------------------------------

strip_tif_extension() {
    local name="$1"
    name="${name%.tif}"
    name="${name%.TIF}"
    name="${name%.tiff}"
    name="${name%.TIFF}"
    printf '%s\n' "$name"
}

extract_parameter() {
    local filename="$1"
    local parameter="$2"

    "$PYTHON" - "$filename" "$parameter" <<'PY'
import os
import re
import sys

name = os.path.basename(sys.argv[1])
parameter = sys.argv[2]
m = re.search(
    rf'(?:^|_){re.escape(parameter)}=([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)',
    name,
)
if m is None:
    raise SystemExit(f"Cannot extract parameter {parameter!r} from filename: {name}")
print(m.group(1))
PY
}

extract_run_number() {
    local run_dir="$1"
    local name="$(basename "$run_dir")"
    if [[ "$name" =~ _run=(.+)$ ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
    else
        die "Cannot extract run number from $run_dir"
    fi
}

prepare_run_output_dirs() {
    local run_dir="$1"
    local d
    log "Preparing output directories: $(basename "$run_dir")"
    for d in \
        "$run_dir/2processed" \
        "$run_dir/3rotated" \
        "$run_dir/4coordinates" \
        "$run_dir/5th2f" \
        "$run_dir/6luminosity"
    do
        rm -rf "$d"
        mkdir -p "$d"
    done
}

# -----------------------------------------------------------------------------
# STEP 1 - LOCAL IMAGE PREPARATION
# -----------------------------------------------------------------------------

cleanup_images() {
    local run_dir="$1"
    local originals_dir="$run_dir/1originals"
    local processed_dir="$run_dir/2processed"
    local input_file filename stem output_file
    local n_diff=0 n_light=0

    log "Local image cleanup (data=diffe is intentionally untouched): $(basename "$run_dir")"

    while IFS= read -r -d '' input_file; do
        filename="$(basename "$input_file")"
        [[ "$filename" == *"data=diffe"* ]] && continue
        if [[ "$filename" != *"data=diff"* && "$filename" != *"data=light"* ]]; then
            continue
        fi

        stem="$(strip_tif_extension "$filename")"
        output_file="$processed_dir/${stem}_processed.tif"
        "$PYTHON" "$PROCESS_IMAGE" \
            --input "$input_file" \
            --process remove_column_bias remove_hot_pixels remove_cold_pixels \
            --output "$output_file"

        if [[ "$filename" == *"data=light"* ]]; then
            ((n_light += 1))
        else
            ((n_diff += 1))
        fi
    done < <(
        find "$originals_dir" -maxdepth 1 -type f \
            \( -iname '*.tif' -o -iname '*.tiff' \) -print0
    )

    (( n_diff > 0 )) || die "No data=diff TIF files found in $originals_dir"
    (( n_light == 1 )) || die "Expected exactly one data=light TIF in $originals_dir; found $n_light"

    if compgen -G "$processed_dir/*data=diffe*" >/dev/null; then
        die "Internal error: a processed data=diffe file exists in $processed_dir"
    fi
}

# -----------------------------------------------------------------------------
# STEP 2 - GEOMETRICAL ALIGNMENT HELPERS
# Used ONLY after local luminosities/systematics have been measured.
# -----------------------------------------------------------------------------

get_unique_light_image() {
    local directory="$1"
    local matches=()
    local f
    while IFS= read -r -d '' f; do matches+=("$f"); done < <(
        find "$directory" -maxdepth 1 -type f \
            \( -iname '*.tif' -o -iname '*.tiff' \) \
            -name '*data=light*_processed.tif' -print0
    )
    (( ${#matches[@]} == 1 )) || die \
        "Expected exactly one processed data=light image in $directory; found ${#matches[@]}."
    printf '%s\n' "${matches[0]}"
}

measure_rotation_angle() {
    local light_image="$1"
    local output angle
    output="$("$PYTHON" "$MEASURE_ROTATION" --input "$light_image")"
    printf '%s\n' "$output" >&2
    angle="$(printf '%s\n' "$output" | "$PYTHON" -c '
import re,sys
text=sys.stdin.read()
m=re.search(r"rotation angle:\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",text)
if m is None: raise SystemExit("Cannot parse rotation angle from measure-rotation.py output.")
print(m.group(1))
')"
    printf '%s\n' "$angle"
}

measure_rotation_clip() {
    local input_image="$1" rotated_image="$2"
    "$PYTHON" - "$input_image" "$rotated_image" <<'PY'
import sys,tifffile
before=tifffile.imread(sys.argv[1]); after=tifffile.imread(sys.argv[2])
ny0,nx0=before.shape[-2:]; ny1,nx1=after.shape[-2:]
print(f"{max(ny0-ny1,0)/2.0:g} {max(nx0-nx1,0)/2.0:g}")
PY
}

rotate_one_image() {
    local input_file="$1" angle="$2" output_file="$3"
    mkdir -p "$(dirname "$output_file")"
    "$PYTHON" "$ROTATE_IMAGE" --input "$input_file" --angle "$angle" --output "$output_file"
}

shift_one_image() {
    local input_file="$1" shifty="$2" shiftx="$3" output_file="$4"
    mkdir -p "$(dirname "$output_file")"
    "$PYTHON" "$SHIFT_IMAGE" --input "$input_file" --shift "$shifty" "$shiftx" --output "$output_file"
}

assert_same_image_shape() {
    local reference="$1" moving="$2"
    "$PYTHON" - "$reference" "$moving" <<'PY'
import sys,tifffile
ref=tifffile.imread(sys.argv[1]); mov=tifffile.imread(sys.argv[2])
if ref.shape != mov.shape:
    raise SystemExit(f"Reference and moving light images differ in shape: {ref.shape} vs {mov.shape}")
PY
}

measure_shift_values() {
    local reference="$1" moving="$2" output
    output="$("$PYTHON" "$MEASURE_SHIFT" --input "$reference" "$moving")"
    printf '%s\n' "$output" >&2
    printf '%s\n' "$output" | "$PYTHON" -c '
import re,sys
text=sys.stdin.read(); number=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
m=re.search(rf"detected subpixel offset\s*\(y,\s*x\)\s*:\s*\[\s*({number})[\s,]+({number})\s*\]",text)
if m is None: raise SystemExit("Cannot parse shift from measure-shift.py output.")
print(m.group(1),m.group(2))
'
}

append_alignment_row() {
    local key="$1" run_dir="$2" angle="$3" clipy="$4" clipx="$5" shifty="$6" shiftx="$7" primary="$8"
    printf '%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$key" "$(basename "$run_dir")" "$angle" "$clipy" "$clipx" \
        "$shifty" "$shiftx" "$primary" >> "$ALIGNMENT_TABLE"
}

# -----------------------------------------------------------------------------
# STEP 3 - LOCAL REFERENCE IMAGE + HOTSPOT COORDINATES
# -----------------------------------------------------------------------------

select_maximum_reference() {
    local processed_dir="$1" varying_parameter="$2"
    "$PYTHON" - "$processed_dir" "$varying_parameter" <<'PY'
from pathlib import Path
import re,sys
d=Path(sys.argv[1]); p=sys.argv[2]
files=[f for f in sorted(d.glob("*data=diff*_processed.tif")) if "data=diffe" not in f.name]
pat=re.compile(rf"(?:^|_){re.escape(p)}=([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)")
e=[]
for f in files:
    m=pat.search(f.name)
    if m: e.append((float(m.group(1)),f))
if not e: raise SystemExit(f"No processed local data=diff images with {p}=... in {d}")
mx=max(v for v,_ in e); matches=[f for v,f in e if v==mx]
if len(matches)!=1: raise SystemExit(f"Expected one maximum-{p} reference; found {len(matches)}")
print(matches[0])
PY
}

rewrite_coordinates_fixed_radius() {
    local source="$1" target="$2" radius="$3"
    "$PYTHON" - "$source" "$target" "$radius" <<'PY'
from pathlib import Path
import math,sys
s=Path(sys.argv[1]); t=Path(sys.argv[2]); r=float(sys.argv[3]); a=math.pi*r*r; n=0
with s.open() as src,t.open('w') as dst:
    for ln,raw in enumerate(src,1):
        if not raw.strip(): continue
        f=[x.strip() for x in raw.split(',')]
        if len(f)<2: raise SystemExit(f"Invalid coordinate line {ln} in {s}")
        x,y=float(f[0]),float(f[1]); dst.write(f"{x:.6f}, {y:.6f}, {a:.6f}, {r:.6f}\n"); n+=1
if n==0: raise SystemExit(f"Coordinate file is empty: {s}")
PY
}

generate_reference_coordinates() {
    local run_dir="$1" key="$2" varying_parameter="$3"
    local coordinates_dir="$run_dir/4coordinates" ref_image reference_coords fixed_r20_coords
    log "Hotspot finding in LOCAL frame: $(basename "$run_dir")" >&2
    ref_image="$(select_maximum_reference "$run_dir/2processed" "$varying_parameter")"
    printf '  local reference: %s\n  selected at maximum %s\n' "$ref_image" "$varying_parameter" >&2
    reference_coords="$coordinates_dir/${key}_detected_reference_coordinates.txt"
    fixed_r20_coords="$coordinates_dir/${key}_detected_R20_coordinates.txt"
    "$PYTHON" "$FIND_DEFECTS" --input "$ref_image" --coordinates_root "$reference_coords" >&2
    [[ -s "$reference_coords" ]] || die "Source finder did not create coordinates: $reference_coords"
    rewrite_coordinates_fixed_radius "$reference_coords" "$fixed_r20_coords" "$R_NOMINAL"
    printf '%s|%s|%s\n' "$reference_coords" "$fixed_r20_coords" "$ref_image"
}

# -----------------------------------------------------------------------------
# STEP 4 - LOCAL TIF -> TH2F
# signal = cleaned local data=diff; error = ORIGINAL local data=diffe
# -----------------------------------------------------------------------------

find_original_diffe() {
    local processed_diff="$1" originals_dir="$2" filename prefix
    filename="$(basename "$processed_diff")"
    prefix="${filename%%data=diff*}"
    local matches=() f
    while IFS= read -r -d '' f; do matches+=("$f"); done < <(
        find "$originals_dir" -maxdepth 1 -type f \
            \( -iname "${prefix}data=diffe*.tif" -o -iname "${prefix}data=diffe*.tiff" \) -print0
    )
    (( ${#matches[@]} == 1 )) || die \
        "Expected exactly one ORIGINAL data=diffe for $filename; found ${#matches[@]} in $originals_dir"
    printf '%s\n' "${matches[0]}"
}

convert_to_th2f() {
    local run_dir="$1" processed_dir="$run_dir/2processed" originals_dir="$run_dir/1originals" root_dir="$run_dir/5th2f"
    local input_file error_file filename stem output_file n_converted=0
    log "Local TIF -> TH2F: cleaned diff + ORIGINAL diffe | $(basename "$run_dir")"
    while IFS= read -r -d '' input_file; do
        [[ "$input_file" == *"data=diffe"* ]] && continue
        error_file="$(find_original_diffe "$input_file" "$originals_dir")"
        filename="$(basename "$input_file")"; stem="$(strip_tif_extension "$filename")"
        output_file="$root_dir/${stem}_th2f.root"
        "$PYTHON" "$TIF2TH2" --input "$input_file" --error "$error_file" --output "$output_file"
        ((n_converted += 1))
    done < <(
        find "$processed_dir" -maxdepth 1 -type f -name '*data=diff*_processed.tif' ! -name '*data=diffe*' -print0
    )
    (( n_converted > 0 )) || die "No local differential images converted in $processed_dir"
}

# -----------------------------------------------------------------------------
# STEP 5 - LOCAL LUMINOSITY AT ONE FIXED RADIUS
# -----------------------------------------------------------------------------

calculate_local_luminosity_radius() {
    local run_dir="$1"
    local key="$2"
    local varying_parameter="$3"
    local source_coords="$4"
    local radius="$5"
    local target_data="$6"

    local run_name
    run_name="$(basename "$run_dir")"

    log "Local luminosity: $run_name | R=$radius px"

    RUN_DIR_ENV="$run_dir" \
    KEY_ENV="$key" \
    VARY_ENV="$varying_parameter" \
    SOURCE_COORDS_ENV="$source_coords" \
    RADIUS_ENV="$radius" \
    TARGET_DATA_ENV="$target_data" \
    DATA_DIR_ENV="$DATA_DIR" \
    SPOT_LUM_DIR_ENV="$SPOT_LUM_DIR" \
    SPOT_LUM_MACRO_ENV="$SPOT_LUM_MACRO" \
    COORD_MATCH_RADIUS_ENV="$COORD_MATCH_RADIUS" \
    VBD_ENV="$VBD_A1" \
    "$PYTHON" <<'PY_LOCAL_RADIUS'
from __future__ import annotations

from pathlib import Path
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

run_dir = Path(os.environ["RUN_DIR_ENV"])
key = os.environ["KEY_ENV"]
vary = os.environ["VARY_ENV"]
source_coords = Path(os.environ["SOURCE_COORDS_ENV"])
radius = float(os.environ["RADIUS_ENV"])
target_data = Path(os.environ["TARGET_DATA_ENV"])
nominal_data = Path(os.environ["DATA_DIR_ENV"])
spot_lum_dir = Path(os.environ["SPOT_LUM_DIR_ENV"])
spot_lum_macro = os.environ["SPOT_LUM_MACRO_ENV"]
coord_tol = float(os.environ["COORD_MATCH_RADIUS_ENV"])
vbd = float(os.environ["VBD_ENV"])
area = math.pi * radius * radius
run_name = run_dir.name
run_number = run_name.split("_run=", 1)[1]


def extract_parameter(filename: str, par: str) -> float:
    m = re.search(
        rf"(?:^|_){re.escape(par)}="
        rf"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
        filename,
    )
    if m is None:
        raise RuntimeError(f"Cannot extract {par}=... from {filename}")
    return float(m.group(1))


def read_local_plan(path: Path) -> pd.DataFrame:
    rows = []
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            fields = [x.strip() for x in line.split(",")]
            if len(fields) < 2:
                raise RuntimeError(f"Invalid coordinate row in {path}: {raw.rstrip()}")
            rows.append({
                "local_spot": len(rows),
                "x": float(fields[0]),
                "y": float(fields[1]),
            })
    if not rows:
        raise RuntimeError(f"Empty coordinate file: {path}")
    return pd.DataFrame(rows)


def ensure_link(link: Path, source: Path):
    if link.is_symlink():
        try:
            if link.resolve() == source.resolve():
                return
        except FileNotFoundError:
            pass
        link.unlink()
    elif link.exists():
        if link.is_dir():
            shutil.rmtree(link)
        else:
            link.unlink()
    link.symlink_to(source.resolve(), target_is_directory=source.is_dir())


# Nominal R=20 lives in DATA itself. R=16/R=24 get light-weight parallel trees.
if target_data.resolve() == nominal_data.resolve():
    target_run = run_dir
    coord_dir = target_run / "4coordinates"
    lum_dir = target_run / "6luminosity"
else:
    target_run = target_data / run_name
    target_run.mkdir(parents=True, exist_ok=True)
    for dirname in ("1originals", "2processed", "3rotated", "5th2f"):
        ensure_link(target_run / dirname, run_dir / dirname)
    coord_dir = target_run / "4coordinates"
    lum_dir = target_run / "6luminosity"
    for d in (coord_dir, lum_dir):
        if d.is_symlink():
            d.unlink()
        elif d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True, exist_ok=True)

coord_dir.mkdir(parents=True, exist_ok=True)
lum_dir.mkdir(parents=True, exist_ok=True)
(target_data / "merged_files" / "local").mkdir(parents=True, exist_ok=True)

# IMPORTANT: these are the coordinates detected in THIS dataset.
plan = read_local_plan(source_coords)
coord_file = coord_dir / f"{key}_detected_R{radius:g}_coordinates.txt"
with coord_file.open("w", encoding="utf-8") as f:
    for row in plan.itertuples(index=False):
        f.write(f"{row.x:.6f}, {row.y:.6f}, {area:.6f}, {radius:.6f}\n")

root_files = sorted(
    p for p in (run_dir / "5th2f").glob("*data=diff*_processed_th2f.root")
    if "data=diffe" not in p.name
)
if not root_files:
    raise RuntimeError(f"No TH2F files found in {run_dir / '5th2f'}")

frames = []
for root_file in root_files:
    t_value = extract_parameter(root_file.name, "T")
    v_value = extract_parameter(root_file.name, "v")
    result_csv = spot_lum_dir / "luminosity_results.csv"
    if result_csv.exists():
        result_csv.unlink()

    macro_call = f'{spot_lum_macro}("{root_file}","{coord_file}")'
    proc = subprocess.run(
        ["root", "-l", "-q", macro_call],
        cwd=spot_lum_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if proc.returncode != 0:
        print(proc.stdout, file=sys.stderr)
        raise RuntimeError(f"ROOT luminosity macro failed for {root_file}")
    if not result_csv.is_file():
        print(proc.stdout, file=sys.stderr)
        raise RuntimeError(f"ROOT did not create luminosity_results.csv for {root_file}")

    result = pd.read_csv(result_csv)
    result_csv.unlink()
    needed = {"x", "y", "luminosity", "error"}
    if not needed.issubset(result.columns):
        raise RuntimeError(
            f"ROOT output for {root_file} is missing columns {sorted(needed-set(result.columns))}"
        )
    for c in needed:
        result[c] = pd.to_numeric(result[c], errors="raise")

    if len(result) != len(plan):
        raise RuntimeError(
            f"ROOT returned {len(result)} luminosity rows for {root_file}; "
            f"{len(plan)} local hotspot coordinates were supplied."
        )

    # Reconnect every ROOT row to the local detected hotspot centre.
    pairs = []
    for ri, r in result.reset_index(drop=True).iterrows():
        for pi, p in plan.iterrows():
            d = math.hypot(float(r.x)-float(p.x), float(r.y)-float(p.y))
            if d <= coord_tol:
                pairs.append((d, ri, pi))
    pairs.sort(key=lambda z: z[0])
    used_r, used_p, mapping = set(), set(), {}
    for d, ri, pi in pairs:
        if ri in used_r or pi in used_p:
            continue
        mapping[ri] = pi
        used_r.add(ri)
        used_p.add(pi)
    if len(mapping) != len(result):
        raise RuntimeError(
            f"Could not reconnect every ROOT result to local coordinates for {root_file}; "
            f"COORD_MATCH_RADIUS={coord_tol} px."
        )

    records = []
    result = result.reset_index(drop=True)
    for ri, r in result.iterrows():
        p = plan.iloc[mapping[ri]]
        records.append({
            "spot": int(p.local_spot),          # local ID until the FINAL stage
            "local_spot": int(p.local_spot),
            "x": float(p.x),
            "y": float(p.y),
            "luminosity": float(r.luminosity),
            "error": float(r.error),
            "T": t_value,
            "v": v_value,
            "v_fin": (v_value - vbd) if key == "T20" else np.nan,
            "dataset_key": key,
            "run_name": run_name,
            "run_number": run_number,
            "detected": True,
            "integration_radius": radius,
            "integration_area": area,
        })

    measured = pd.DataFrame(records)
    measured = measured.sort_values(["local_spot", vary], kind="stable").reset_index(drop=True)
    point_output = lum_dir / f"local_R{radius:g}_luminosity_{root_file.stem}.csv"
    measured.to_csv(point_output, index=False)
    frames.append(measured)

run_master = pd.concat(frames, ignore_index=True, sort=False)
run_master = run_master.sort_values(["local_spot", vary], kind="stable").reset_index(drop=True)

run_output = lum_dir / f"{run_name}_local_R{radius:g}.csv"
merged_local = target_data / "merged_files" / "local" / f"{run_name}_local_R{radius:g}.csv"
run_master.to_csv(run_output, index=False)
run_master.to_csv(merged_local, index=False)

print(f"Created local R={radius:g} run file: {run_output}")
print(f"Local hotspot count: {run_master['local_spot'].nunique()}")
PY_LOCAL_RADIUS
}

# -----------------------------------------------------------------------------
# PROCESS EACH DATASET LOCALLY THROUGH SOURCE FINDING + TH2F
# -----------------------------------------------------------------------------

process_local_run() {
    local run_dir="$1" key="$2" varying_parameter="$3"
    printf '\n============================================================\n'
    printf 'LOCAL MEASUREMENT PREPARATION: %s\n' "$(basename "$run_dir")"
    printf 'HOTSPOT FINDING: %s at maximum %s\n' "$key" "$varying_parameter"
    printf '============================================================\n'
    prepare_run_output_dirs "$run_dir"
    cleanup_images "$run_dir"
    local coord_triplet source_coords
    coord_triplet="$(generate_reference_coordinates "$run_dir" "$key" "$varying_parameter")"
    source_coords="${coord_triplet%%|*}"
    convert_to_th2f "$run_dir"
    LAST_SOURCE_COORDS="$source_coords"
}

# -----------------------------------------------------------------------------
# STEP 6 - SYSTEMATIC UNCERTAINTY WHILE IDs ARE STILL LOCAL
# -----------------------------------------------------------------------------

calculate_local_systematic_uncertainty() {
    log "Local luminosity systematic uncertainty from R=16/20/24"

    MANIFEST_ENV="$MANIFEST" \
    DATA20_ENV="$DATA_DIR" \
    DATA16_ENV="$R16_DATA_DIR" \
    DATA24_ENV="$R24_DATA_DIR" \
    "$PYTHON" <<'PY_LOCAL_SYS'
from pathlib import Path
import os
import shutil
import numpy as np
import pandas as pd

manifest = pd.read_csv(os.environ["MANIFEST_ENV"])
data = {
    20: Path(os.environ["DATA20_ENV"]),
    16: Path(os.environ["DATA16_ENV"]),
    24: Path(os.environ["DATA24_ENV"]),
}

all_diag = []
key_cols = ["dataset_key", "run_number", "local_spot", "T", "v"]

for m in manifest.itertuples(index=False):
    paths = {
        r: data[r] / "merged_files" / "local" / f"{m.run_name}_local_R{r}.csv"
        for r in (16, 20, 24)
    }
    for p in paths.values():
        if not p.is_file():
            raise RuntimeError(f"Required local-radius file not found: {p}")

    tables = {r: pd.read_csv(p) for r, p in paths.items()}
    for r, df in tables.items():
        missing = set(key_cols + ["x", "y", "luminosity", "error", "integration_radius"]) - set(df.columns)
        if missing:
            raise RuntimeError(f"R={r} local file is missing columns {sorted(missing)}: {paths[r]}")
        if df.duplicated(key_cols).any():
            raise RuntimeError(f"R={r} local file has duplicate measurement keys: {paths[r]}")
        if not np.allclose(pd.to_numeric(df["integration_radius"]), float(r), atol=1e-6, rtol=0):
            raise RuntimeError(f"Unexpected integration radius in {paths[r]}")

    idx = {r: tables[r].set_index(key_cols) for r in (16, 20, 24)}
    if set(idx[20].index) != set(idx[16].index) or set(idx[20].index) != set(idx[24].index):
        raise RuntimeError(
            f"R=16/20/24 do not contain the same LOCAL measurements for {m.run_name}"
        )

    # The same locally detected centre MUST be used for all three radii.
    for k in idx[20].index:
        p16, p20, p24 = idx[16].loc[k], idx[20].loc[k], idx[24].loc[k]
        if max(
            abs(float(p16.x)-float(p20.x)),
            abs(float(p24.x)-float(p20.x)),
            abs(float(p16.y)-float(p20.y)),
            abs(float(p24.y)-float(p20.y)),
        ) > 1e-3:
            raise RuntimeError(f"Local coordinate mismatch across radii for {m.run_name}, key={k}")

    L16 = idx[16]["luminosity"].astype(float)
    L20 = idx[20]["luminosity"].astype(float)
    L24 = idx[24]["luminosity"].astype(float)
    d1 = (L24-L20).abs()
    d2 = (L20-L16).abs()
    d = pd.concat([d1, d2], axis=1).max(axis=1)

    unc = pd.DataFrame({"deltaL1": d1, "deltaL2": d2, "deltaL": d}).reset_index()
    nominal = tables[20].drop(columns=[c for c in ("deltaL1", "deltaL2", "deltaL") if c in tables[20].columns])
    nominal = nominal.merge(unc, on=key_cols, how="left", validate="one_to_one")
    if nominal[["deltaL1", "deltaL2", "deltaL"]].isna().any().any():
        raise RuntimeError(f"Failed to attach every local systematic uncertainty for {m.run_name}")

    # Update both nominal local copies: run/6luminosity and merged_files/local.
    targets = [
        Path(m.run_dir) / "6luminosity" / f"{m.run_name}_local_R20.csv",
        paths[20],
    ]
    for target in targets:
        backup = Path(str(target) + ".pre_systematic_uncertainty.bak")
        if target.is_file() and not backup.exists():
            shutil.copy2(target, backup)
        nominal.to_csv(target, index=False)

    diag = idx[20].reset_index()[key_cols + ["x", "y", "luminosity", "error"]].rename(
        columns={"luminosity": "L20", "error": "stat_error_R20"}
    )
    diag = diag.merge(
        idx[16].reset_index()[key_cols + ["luminosity"]].rename(columns={"luminosity":"L16"}),
        on=key_cols, validate="one_to_one"
    ).merge(
        idx[24].reset_index()[key_cols + ["luminosity"]].rename(columns={"luminosity":"L24"}),
        on=key_cols, validate="one_to_one"
    ).merge(unc, on=key_cols, validate="one_to_one")
    diag["deltaL_over_absL20"] = np.where(diag["L20"].abs()>0, diag["deltaL"]/diag["L20"].abs(), np.nan)
    all_diag.append(diag)

local_diag = pd.concat(all_diag, ignore_index=True, sort=False)
out = data[20] / "merged_files" / "A1_local_systematic_uncertainty_diagnostics.csv"
local_diag.to_csv(out, index=False)
print(f"Created local systematic diagnostic table: {out}")
print("Systematic uncertainties were computed BEFORE global-ID assignment.")
PY_LOCAL_SYS
}

# -----------------------------------------------------------------------------
# STEP 7 - BUILD COMMON MATCHING COORDINATES (AFTER ALL MEASUREMENTS)
# -----------------------------------------------------------------------------

transform_coordinates_for_matching() {
    local run_dir="$1" key="$2" source_coords="$3" local_reference="$4"
    local angle="$5" shifty="$6" shiftx="$7" is_primary="$8"
    local output_csv="$run_dir/4coordinates/${key}_matching_aligned_coordinates.csv"

    RUN_DIR_ENV="$run_dir" KEY_ENV="$key" SOURCE_COORDS_ENV="$source_coords" \
    LOCAL_REFERENCE_ENV="$local_reference" ANGLE_ENV="$angle" \
    SHIFTY_ENV="$shifty" SHIFTX_ENV="$shiftx" PRIMARY_ENV="$is_primary" \
    OUTPUT_ENV="$output_csv" PYTHON_ENV="$PYTHON" ROTATE_ENV="$ROTATE_IMAGE" SHIFT_ENV="$SHIFT_IMAGE" \
    "$PYTHON" <<'PY_MATCH_COORDS'
from pathlib import Path
import math,os,shutil,subprocess
import numpy as np
import pandas as pd
import tifffile

run_dir=Path(os.environ['RUN_DIR_ENV'])
coords_path=Path(os.environ['SOURCE_COORDS_ENV'])
ref_path=Path(os.environ['LOCAL_REFERENCE_ENV'])
angle=float(os.environ['ANGLE_ENV']); shifty=float(os.environ['SHIFTY_ENV']); shiftx=float(os.environ['SHIFTX_ENV'])
is_primary=os.environ['PRIMARY_ENV'].lower() in ('1','true','yes','y')
out_path=Path(os.environ['OUTPUT_ENV'])
py=os.environ['PYTHON_ENV']; rotate=os.environ['ROTATE_ENV']; shift=os.environ['SHIFT_ENV']

coords=[]
with coords_path.open() as f:
    for raw in f:
        if not raw.strip(): continue
        q=[x.strip() for x in raw.split(',')]
        coords.append((len(coords),float(q[0]),float(q[1])))
if not coords: raise RuntimeError(f'No coordinates in {coords_path}')

ref=tifffile.imread(ref_path)
h,w=ref.shape[-2:]
tmp=run_dir/'.matching_coordinate_transform_tmp'
if tmp.exists(): shutil.rmtree(tmp)
tmp.mkdir(parents=True)
rows=[]
try:
    for local_spot,x_root,y_root in coords:
        x_img=x_root; y_img=h-y_root
        marker=np.zeros((h,w),dtype=np.float32)
        r0=6
        x0=max(0,int(math.floor(x_img))-r0); x1=min(w,int(math.floor(x_img))+r0+2)
        y0=max(0,int(math.floor(y_img))-r0); y1=min(h,int(math.floor(y_img))+r0+2)
        sy,sx=np.mgrid[y0:y1,x0:x1]
        marker[y0:y1,x0:x1]=np.exp(-((sx-x_img)**2+(sy-y_img)**2)/(2*1.2**2)).astype(np.float32)

        inp=tmp/f'marker_{local_spot:04d}_input.tif'
        rot=tmp/f'marker_{local_spot:04d}_rotated.tif'
        final=rot if is_primary else tmp/f'marker_{local_spot:04d}_aligned.tif'
        tifffile.imwrite(inp,marker)
        subprocess.run([py,rotate,'--input',str(inp),'--angle',str(angle),'--output',str(rot)],check=True,
                       stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True)
        if not is_primary:
            subprocess.run([py,shift,'--input',str(rot),'--shift',str(shifty),str(shiftx),'--output',str(final)],check=True,
                           stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True)

        arr=np.asarray(tifffile.imread(final),dtype=float)
        wt=np.abs(arr); total=wt.sum()
        if not np.isfinite(total) or total<=0:
            raise RuntimeError(f'Transformed marker vanished for local_spot={local_spot}')
        yy,xx=np.indices(wt.shape[-2:])
        x_match=float((wt*xx).sum()/total)
        y_img_match=float((wt*yy).sum()/total)
        y_match=float(wt.shape[-2]-y_img_match)
        rows.append({'local_spot':local_spot,'x_local':x_root,'y_local':y_root,
                     'x_match':x_match,'y_match':y_match})
finally:
    shutil.rmtree(tmp,ignore_errors=True)

out=pd.DataFrame(rows).sort_values('local_spot',kind='stable')
out_path.parent.mkdir(parents=True,exist_ok=True)
out.to_csv(out_path,index=False)
print(f'Created matching-coordinate table: {out_path}')
PY_MATCH_COORDS
    printf '%s\n' "$output_csv"
}

prepare_matching_geometry() {
    local coords_v9="$1" coords_v7="$2" coords_t20="$3"
    log "Geometrical alignment for GLOBAL-ID MATCHING ONLY"
    printf '%s\n' \
        "dataset_key,run_name,angle,clipy,clipx,shifty,shiftx,is_primary_reference" \
        > "$ALIGNMENT_TABLE"

    local run_dir key vary source_coords local_ref processed_light angle
    local rotated_light prealign_light clipy clipx shifty shiftx

    run_dir="$RUN_V9"; key="v9"; vary="T"; source_coords="$coords_v9"
    local_ref="$(select_maximum_reference "$run_dir/2processed" "$vary")"
    processed_light="$(get_unique_light_image "$run_dir/2processed")"
    assert_same_image_shape "$local_ref" "$processed_light"
    angle="$(measure_rotation_angle "$processed_light")"
    rotated_light="$run_dir/3rotated/${key}_light_aligned.tif"
    rotate_one_image "$processed_light" "$angle" "$rotated_light"
    BASE_LIGHT_REFERENCE="$rotated_light"
    read -r clipy clipx < <(measure_rotation_clip "$processed_light" "$rotated_light")
    shifty=0; shiftx=0
    append_alignment_row "$key" "$run_dir" "$angle" "$clipy" "$clipx" "$shifty" "$shiftx" yes
    transform_coordinates_for_matching "$run_dir" "$key" "$source_coords" "$local_ref" \
        "$angle" "$shifty" "$shiftx" yes >/dev/null

    local spec
    for spec in "v7|$RUN_V7|T|$coords_v7" "T20|$RUN_T20|v|$coords_t20"; do
        IFS='|' read -r key run_dir vary source_coords <<< "$spec"
        local_ref="$(select_maximum_reference "$run_dir/2processed" "$vary")"
        processed_light="$(get_unique_light_image "$run_dir/2processed")"
        assert_same_image_shape "$local_ref" "$processed_light"
        angle="$(measure_rotation_angle "$processed_light")"
        prealign_light="$run_dir/3rotated/${key}_light_rotated_pre_shift.tif"
        rotated_light="$run_dir/3rotated/${key}_light_aligned.tif"
        rotate_one_image "$processed_light" "$angle" "$prealign_light"
        read -r clipy clipx < <(measure_rotation_clip "$processed_light" "$prealign_light")
        assert_same_image_shape "$BASE_LIGHT_REFERENCE" "$prealign_light"
        read -r shifty shiftx < <(measure_shift_values "$BASE_LIGHT_REFERENCE" "$prealign_light")
        shift_one_image "$prealign_light" "$shifty" "$shiftx" "$rotated_light"
        rm -f "$prealign_light"
        append_alignment_row "$key" "$run_dir" "$angle" "$clipy" "$clipx" "$shifty" "$shiftx" no
        transform_coordinates_for_matching "$run_dir" "$key" "$source_coords" "$local_ref" \
            "$angle" "$shifty" "$shiftx" no >/dev/null
    done
}

# -----------------------------------------------------------------------------
# STEP 8 - FINAL GLOBAL-ID ASSIGNMENT (NO FORCED MEASUREMENTS)
# -----------------------------------------------------------------------------


assign_global_ids_after_measurement() {
    log "Final global-ID assignment from locally detected hotspots"

    SENSOR_ENV="$SENSOR" \
    MANIFEST_ENV="$MANIFEST" \
    DATA20_ENV="$DATA_DIR" \
    DATA16_ENV="$R16_DATA_DIR" \
    DATA24_ENV="$R24_DATA_DIR" \
    MATCH_RADIUS_ENV="$MATCH_RADIUS" \
    "$PYTHON" <<'PY_GLOBAL_FINAL'
from __future__ import annotations

from pathlib import Path
import math
import os
import numpy as np
import pandas as pd

sensor = os.environ["SENSOR_ENV"]
manifest = pd.read_csv(os.environ["MANIFEST_ENV"])
data_dirs = {
    20: Path(os.environ["DATA20_ENV"]),
    16: Path(os.environ["DATA16_ENV"]),
    24: Path(os.environ["DATA24_ENV"]),
}
match_radius = float(os.environ["MATCH_RADIUS_ENV"])

if list(manifest["dataset_key"]) != ["v9", "v7", "T20"]:
    raise RuntimeError("Manifest order must be v9, v7, T20 so v9 seeds the global IDs.")

# Local x,y are the physical measurement coordinates. x_match,y_match are
# auxiliary coordinates in the common aligned frame and are used only for ID matching.
datasets = []
for order, m in enumerate(manifest.itertuples(index=False)):
    local_path = data_dirs[20] / "merged_files" / "local" / f"{m.run_name}_local_R20.csv"
    local_df = pd.read_csv(local_path)
    local_coords = (local_df[["local_spot","x","y"]].drop_duplicates()
                    .sort_values("local_spot",kind="stable").reset_index(drop=True))
    match_path = Path(m.run_dir) / "4coordinates" / f"{m.dataset_key}_matching_aligned_coordinates.csv"
    if not match_path.is_file():
        raise RuntimeError(f"Matching-coordinate file not found: {match_path}")
    match_df = pd.read_csv(match_path)
    required={"local_spot","x_local","y_local","x_match","y_match"}
    if not required.issubset(match_df.columns):
        raise RuntimeError(f"Missing matching-coordinate columns in {match_path}")
    coords=local_coords.merge(match_df[["local_spot","x_match","y_match"]],
                              on="local_spot",how="left",validate="one_to_one")
    if coords[["x_match","y_match"]].isna().any().any():
        raise RuntimeError(f"Missing transformed coordinates in {match_path}")
    datasets.append({
        "dataset_key":str(m.dataset_key),"run_name":str(m.run_name),
        "run_number":str(m.run_number),"run_dir":Path(m.run_dir),
        "varying_parameter":str(m.varying_parameter),"order":order,"coords":coords,
    })

def greedy_match(source: pd.DataFrame, catalog: pd.DataFrame, radius: float):
    pairs = []
    counts = {int(r.local_spot): 0 for r in source.itertuples(index=False)}
    for s in source.itertuples(index=False):
        for t in catalog.itertuples(index=False):
            dist = math.hypot(float(s.x_match)-float(t.x_match_ref), float(s.y_match)-float(t.y_match_ref))
            if dist <= radius:
                pairs.append((dist, int(s.local_spot), int(t.spot)))
                counts[int(s.local_spot)] += 1
    pairs.sort(key=lambda z: z[0])
    used_local, used_global = set(), set()
    mapping, distances = {}, {}
    for dist, local, gid in pairs:
        if local in used_local or gid in used_global:
            continue
        mapping[local] = gid
        distances[local] = dist
        used_local.add(local)
        used_global.add(gid)
    return mapping, distances, counts


# v9 is only the ID seed. Its coordinates are NOT propagated into other datasets.
primary = datasets[0]
catalog_rows = []
primary_mapping = {}
for gid, row in enumerate(primary["coords"].itertuples(index=False)):
    local = int(row.local_spot)
    primary_mapping[local] = gid
    catalog_rows.append({
        "spot": gid,
        "x_match_ref": float(row.x_match),
        "y_match_ref": float(row.y_match),
        "x_local_seed": float(row.x),
        "y_local_seed": float(row.y),
        "first_seen_dataset": primary["dataset_key"],
        "first_seen_run": primary["run_number"],
    })

catalog = pd.DataFrame(catalog_rows)
next_gid = len(catalog)
mapping_rows = []

for d in datasets:
    if d is primary:
        mapping = dict(primary_mapping)
        distances = {local: 0.0 for local in mapping}
        counts = {local: 1 for local in mapping}
        types = {local: "primary_reference" for local in mapping}
    else:
        mapping, distances, counts = greedy_match(d["coords"], catalog, match_radius)
        types = {local: "matched" for local in mapping}

        # A hotspot absent from previous datasets simply receives a new global ID.
        # It is NOT force-measured in those datasets.
        for row in d["coords"].itertuples(index=False):
            local = int(row.local_spot)
            if local in mapping:
                continue
            gid = next_gid
            next_gid += 1
            mapping[local] = gid
            distances[local] = np.nan
            counts.setdefault(local, 0)
            types[local] = "new"
            catalog = pd.concat([
                catalog,
                pd.DataFrame([{
                    "spot": gid,
                    "x_match_ref": float(row.x_match),
                    "y_match_ref": float(row.y_match),
                    "x_local_seed": float(row.x),
                    "y_local_seed": float(row.y),
                    "first_seen_dataset": d["dataset_key"],
                    "first_seen_run": d["run_number"],
                }])
            ], ignore_index=True)

    d["mapping"] = mapping
    for row in d["coords"].itertuples(index=False):
        local = int(row.local_spot)
        mapping_rows.append({
            "dataset_key": d["dataset_key"],
            "run_name": d["run_name"],
            "run_number": d["run_number"],
            "local_spot": local,
            "spot": int(mapping[local]),
            "x_detected": float(row.x),
            "y_detected": float(row.y),
            "x_match": float(row.x_match),
            "y_match": float(row.y_match),
            "match_distance": distances.get(local, np.nan),
            "match_type": types[local],
            "candidate_count": int(counts.get(local, 0)),
        })

catalog["spot"] = catalog["spot"].astype(int)
catalog = catalog.sort_values("spot", kind="stable").reset_index(drop=True)
mapping_df = pd.DataFrame(mapping_rows).sort_values(["spot", "dataset_key"], kind="stable")

merged20 = data_dirs[20] / "merged_files"
merged20.mkdir(parents=True, exist_ok=True)
catalog_path = merged20 / f"{sensor}_global_catalog.csv"
mapping_path = merged20 / f"{sensor}_global_mapping.csv"
catalog.to_csv(catalog_path, index=False)
mapping_df.to_csv(mapping_path, index=False)
print(f"Created: {catalog_path}")
print(f"Created: {mapping_path}")

# Apply the SAME local->global mapping to each radius. No rows are added.
for radius in (20, 16, 24):
    data_dir = data_dirs[radius]
    merged = data_dir / "merged_files"
    merged.mkdir(parents=True, exist_ok=True)
    all_runs = []

    # Keep catalog/mapping copies in each radius tree for diagnostics.
    if radius != 20:
        catalog.to_csv(merged / f"{sensor}_global_catalog.csv", index=False)
        mapping_df.to_csv(merged / f"{sensor}_global_mapping.csv", index=False)

    for d in datasets:
        local_path = merged / "local" / f"{d['run_name']}_local_R{radius}.csv"
        df = pd.read_csv(local_path)
        map_series = pd.Series(d["mapping"], name="spot")
        df["spot"] = df["local_spot"].map(map_series)
        if df["spot"].isna().any():
            raise RuntimeError(f"Unmapped local hotspots remain in {local_path}")
        df["spot"] = df["spot"].astype(int)
        df["detected"] = True

        md = mapping_df[mapping_df["dataset_key"] == d["dataset_key"]][
            ["local_spot", "x_match", "y_match", "match_distance", "match_type"]
        ]
        df = df.drop(columns=[c for c in ("x_match", "y_match", "match_distance", "match_type") if c in df.columns])
        df = df.merge(md, on="local_spot", how="left", validate="many_to_one")

        vary = d["varying_parameter"]
        df = df.sort_values(["spot", vary], kind="stable").reset_index(drop=True)

        target_run = d["run_dir"] if radius == 20 else data_dir / d["run_name"]
        lum_dir = target_run / "6luminosity"
        lum_dir.mkdir(parents=True, exist_ok=True)
        run_output = lum_dir / f"{d['run_name']}_global_ID.csv"
        merged_output = merged / f"{d['run_name']}_global_ID.csv"
        df.to_csv(run_output, index=False)
        df.to_csv(merged_output, index=False)
        print(f"Created R={radius} global-ID file: {merged_output}")

        df["__dataset_order"] = d["order"]
        all_runs.append(df)

    master = pd.concat(all_runs, ignore_index=True, sort=False)
    master = master.sort_values(["spot", "__dataset_order", "T", "v"], kind="stable").reset_index(drop=True)
    master = master.drop(columns=["__dataset_order"])
    master_path = merged / f"{sensor}_all_runs_global_ID.csv"
    master.to_csv(master_path, index=False)
    print(f"Created R={radius} master: {master_path}")

# Final diagnostics expressed with global IDs, based only on actually measured rows.
m20 = pd.read_csv(data_dirs[20] / "merged_files" / f"{sensor}_all_runs_global_ID.csv")
m16 = pd.read_csv(data_dirs[16] / "merged_files" / f"{sensor}_all_runs_global_ID.csv")
m24 = pd.read_csv(data_dirs[24] / "merged_files" / f"{sensor}_all_runs_global_ID.csv")
keys = ["dataset_key", "run_number", "spot", "T", "v"]
for label, df in (("R20",m20),("R16",m16),("R24",m24)):
    if df.duplicated(keys).any():
        raise RuntimeError(f"Duplicate final measurement keys in {label}")
idx16, idx20, idx24 = (df.set_index(keys) for df in (m16,m20,m24))
if set(idx20.index) != set(idx16.index) or set(idx20.index) != set(idx24.index):
    raise RuntimeError("Final radius masters differ in their actually measured rows")

diag = idx20.reset_index()[keys + ["local_spot", "x", "y", "luminosity", "error", "deltaL1", "deltaL2", "deltaL"]].rename(
    columns={"luminosity":"L20", "error":"stat_error_R20"}
)
diag = diag.merge(
    idx16.reset_index()[keys+["luminosity"]].rename(columns={"luminosity":"L16"}),
    on=keys, validate="one_to_one"
).merge(
    idx24.reset_index()[keys+["luminosity"]].rename(columns={"luminosity":"L24"}),
    on=keys, validate="one_to_one"
)
diag["deltaL_over_absL20"] = np.where(diag["L20"].abs()>0, diag["deltaL"]/diag["L20"].abs(), np.nan)
diag_path = data_dirs[20] / "merged_files" / f"{sensor}_systematic_uncertainty_diagnostics.csv"
diag.to_csv(diag_path, index=False)
print(f"Created final systematic diagnostic table: {diag_path}")

print("")
print("GLOBAL-ID ASSIGNMENT COMPLETED")
print(f"Union of detected global hotspots: {len(catalog)}")
for d in datasets:
    print(f"  {d['dataset_key']}: {len(d['mapping'])} detected hotspots")
print("No forced measurements were created.")
PY_GLOBAL_FINAL
}

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

main() {
    check_requirements
    mkdir -p "$MERGED_DIR" "$R16_DATA_DIR/merged_files" "$R24_DATA_DIR/merged_files"

    printf '%s\n' \
        "dataset_key,run_name,run_number,varying_parameter,run_dir,detected_coordinates" \
        > "$MANIFEST"

    printf '%s\n' \
        "dataset_key,run_name,angle,clipy,clipx,shifty,shiftx,is_primary_reference" \
        > "$ALIGNMENT_TABLE"

    printf '\n============================================================\n'
    printf 'REDUCED GOLD-SAMPLE EMMI PIPELINE\n'
    printf 'BASE_DIR:          %s\n' "$BASE_DIR"
    printf 'EMMI_DIR:          %s\n' "$EMMI_DIR"
    printf 'DATA_DIR:          %s\n' "$DATA_DIR"
    printf 'SOURCE FINDING:    max T for v=9/v=7; max Vover for T=20\n'
    printf 'LOCAL RADII:       R=20,16,24 use each dataset own detected centres\n'
    printf 'STAT ERROR MAP:    ORIGINAL data=diffe, never cleaned/rotated/shifted\n'
    printf 'GLOBAL IDS:        assigned only after local luminosity + systematics\n'
    printf 'FORCED MEASURE:    disabled\n'
    printf '============================================================\n'

    # Measure all datasets independently in their own native frame.
    local coords_v9 coords_v7 coords_t20

    process_local_run "$RUN_V9" "v9" "T"
    coords_v9="$LAST_SOURCE_COORDS"

    process_local_run "$RUN_V7" "v7" "T"
    coords_v7="$LAST_SOURCE_COORDS"

    process_local_run "$RUN_T20" "T20" "v"
    coords_t20="$LAST_SOURCE_COORDS"

    printf '%s,%s,%s,%s,%s,%s\n' \
        "v9" "$(basename "$RUN_V9")" "$(extract_run_number "$RUN_V9")" "T" "$RUN_V9" "$coords_v9" \
        >> "$MANIFEST"
    printf '%s,%s,%s,%s,%s,%s\n' \
        "v7" "$(basename "$RUN_V7")" "$(extract_run_number "$RUN_V7")" "T" "$RUN_V7" "$coords_v7" \
        >> "$MANIFEST"
    printf '%s,%s,%s,%s,%s,%s\n' \
        "T20" "$(basename "$RUN_T20")" "$(extract_run_number "$RUN_T20")" "v" "$RUN_T20" "$coords_t20" \
        >> "$MANIFEST"

    # ---------------------------------------------------------
    # Luminosities FIRST: local coordinates and local IDs only.
    # ---------------------------------------------------------
    calculate_local_luminosity_radius "$RUN_V9"  "v9"  "T" "$coords_v9"  "$R_NOMINAL" "$DATA_DIR"
    calculate_local_luminosity_radius "$RUN_V7"  "v7"  "T" "$coords_v7"  "$R_NOMINAL" "$DATA_DIR"
    calculate_local_luminosity_radius "$RUN_T20" "T20" "v" "$coords_t20" "$R_NOMINAL" "$DATA_DIR"

    calculate_local_luminosity_radius "$RUN_V9"  "v9"  "T" "$coords_v9"  "$R_LOW" "$R16_DATA_DIR"
    calculate_local_luminosity_radius "$RUN_V7"  "v7"  "T" "$coords_v7"  "$R_LOW" "$R16_DATA_DIR"
    calculate_local_luminosity_radius "$RUN_T20" "T20" "v" "$coords_t20" "$R_LOW" "$R16_DATA_DIR"

    calculate_local_luminosity_radius "$RUN_V9"  "v9"  "T" "$coords_v9"  "$R_HIGH" "$R24_DATA_DIR"
    calculate_local_luminosity_radius "$RUN_V7"  "v7"  "T" "$coords_v7"  "$R_HIGH" "$R24_DATA_DIR"
    calculate_local_luminosity_radius "$RUN_T20" "T20" "v" "$coords_t20" "$R_HIGH" "$R24_DATA_DIR"

    # Systematic uncertainty is still computed with LOCAL IDs.
    calculate_local_systematic_uncertainty

    # Only now build a common matching frame from the light images.
    prepare_matching_geometry "$coords_v9" "$coords_v7" "$coords_t20"

    # Final bookkeeping only: no luminosity is recomputed and no absent spot is forced.
    assign_global_ids_after_measurement

    printf '\n============================================================\n'
    printf 'PIPELINE COMPLETED SUCCESSFULLY\n'
    printf 'Nominal master:          %s\n' "$MERGED_DIR/${SENSOR}_all_runs_global_ID.csv"
    printf 'R=16 master:             %s\n' "$R16_DATA_DIR/merged_files/${SENSOR}_all_runs_global_ID.csv"
    printf 'R=24 master:             %s\n' "$R24_DATA_DIR/merged_files/${SENSOR}_all_runs_global_ID.csv"
    printf 'Systematic diagnostics:  %s\n' "$MERGED_DIR/${SENSOR}_systematic_uncertainty_diagnostics.csv"
    printf 'Alignment table:         %s\n' "$ALIGNMENT_TABLE"
    printf 'Global catalog:          %s\n' "$MERGED_DIR/${SENSOR}_global_catalog.csv"
    printf 'Forced measurements:     NONE\n'
    printf '============================================================\n'
}

main "$@"
