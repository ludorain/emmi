#!/usr/bin/env bash

# ============================================================
# EMMI - RECALCULATE LUMINOSITY WITH AN ARBITRARY FIXED RADIUS
#
# This is a SECOND pipeline. It does NOT repeat image cleanup,
# rotation, alignment, source finding, or TIF -> TH2F conversion.
# It reuses the validated products of:
#
#   emmi/DATA_irradiated_isolated_changeR
#
# and writes the radius-specific analysis to:
#
#   emmi/DATA_irradiated_isolated_R=<R>
#
# Usage:
#   ./recalculate_isolated_radius.sh A1 T 8
#   ./recalculate_isolated_radius.sh A1 v 10.5
#   ./recalculate_isolated_radius.sh B2 T 12
#
# IMPORTANT GLOBAL-ID POLICY
# --------------------------
# Global hotspot IDs are NEVER reassigned in this pipeline.
# They are inherited from the source pipeline's validated
# geometry_plan/global_mapping. Therefore the same physical hotspot
# has exactly the same global ID for every tested integration radius.
# ============================================================

set -Eeuo pipefail
shopt -s nullglob

# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Authoritative source dataset. This directory must already have been
# processed with the full pipeline (including PHASE 3).
SOURCE_DATA_DIR="${SOURCE_DATA_DIR:-$BASE_DIR/DATA_irradiated_isolated_changeR}"

SPOT_LUM_DIR="$BASE_DIR/spot_luminosity"
SPOT_LUM_MACRO="spot_luminosity_sum_irradiated.C"
PYTHON="${PYTHON:-python3}"

# This tolerance is used ONLY to reconnect ROOT output rows to the exact
# coordinate plan supplied to ROOT. It does NOT assign global hotspot IDs.
COORD_MATCH_RADIUS="${COORD_MATCH_RADIUS:-1.0}"

# Breakdown voltages, used only to calculate v_fin when T is fixed.
VBD_A1="${VBD_A1:-51.3}"
VBD_A2="${VBD_A2:-}"
VBD_B1="${VBD_B1:-50.9}"
VBD_B2="${VBD_B2:-}"

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
    printf '\nERROR: radius-recalculation pipeline stopped at line %s (exit code %s).\n' \
        "$line_no" "$exit_code" >&2
    exit "$exit_code"
}

trap 'on_error $LINENO' ERR

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------

if [[ $# -ne 3 ]]; then
    cat >&2 <<USAGE
Usage:
  $0 SENSOR CONSTANT RADIUS

SENSOR:
  A1 | A2 | B1 | B2

CONSTANT:
  T | v

RADIUS:
  positive integration radius in pixels

Examples:
  $0 A1 T 8
  $0 A1 v 10.5
  $0 B2 T 12
USAGE
    exit 1
fi

SENSOR="$1"
CONSTANT="$2"
RADIUS_LABEL="$3"

case "$SENSOR" in
    A1|A2|B1|B2) ;;
    *) die "Invalid sensor '$SENSOR'. Allowed values: A1 A2 B1 B2." ;;
esac

case "$CONSTANT" in
    T|v) ;;
    *) die "Invalid constant '$CONSTANT'. Allowed values: T or v." ;;
esac

# Accept ordinary positive decimal/scientific notation, but reject paths,
# expressions, commas, etc. The original textual value is preserved in the
# output-directory name.
if ! [[ "$RADIUS_LABEL" =~ ^[+]?(\.?[0-9]+|[0-9]+\.?[0-9]*)([eE][-+]?[0-9]+)?$ ]]; then
    die "Invalid radius '$RADIUS_LABEL'. Give one positive numeric value in pixels."
fi

RADIUS_VALUE="$($PYTHON - "$RADIUS_LABEL" <<'PY'
import math
import sys
r = float(sys.argv[1])
if not math.isfinite(r) or r <= 0:
    raise SystemExit(1)
print(repr(r))
PY
)" || die "Radius must be finite and > 0."

TARGET_DATA_DIR="$BASE_DIR/DATA_irradiated_isolated_R=$RADIUS_LABEL"
TARGET_MERGED_ROOT="$TARGET_DATA_DIR/merged_files"

# ------------------------------------------------------------
# REQUIREMENTS
# ------------------------------------------------------------

command -v "$PYTHON" >/dev/null 2>&1 \
    || die "Python executable not found: $PYTHON"
command -v root >/dev/null 2>&1 \
    || die "ROOT executable 'root' not found in PATH."

[[ -d "$SOURCE_DATA_DIR" ]] \
    || die "Source DATA directory not found: $SOURCE_DATA_DIR"
[[ -f "$SPOT_LUM_DIR/$SPOT_LUM_MACRO" ]] \
    || die "ROOT luminosity macro not found: $SPOT_LUM_DIR/$SPOT_LUM_MACRO"

get_vbd() {
    case "$SENSOR" in
        A1) printf '%s\n' "$VBD_A1" ;;
        A2) printf '%s\n' "$VBD_A2" ;;
        B1) printf '%s\n' "$VBD_B1" ;;
        B2) printf '%s\n' "$VBD_B2" ;;
    esac
}

VBD="$(get_vbd)"

mkdir -p "$TARGET_DATA_DIR" "$TARGET_MERGED_ROOT"

printf '\n'
printf '============================================================\n'
printf 'EMMI FIXED-RADIUS RECALCULATION PIPELINE\n'
printf 'Source DATA:       %s\n' "$SOURCE_DATA_DIR"
printf 'Target DATA:       %s\n' "$TARGET_DATA_DIR"
printf 'Sensor:            %s\n' "$SENSOR"
printf 'Constant:          %s\n' "$CONSTANT"
printf 'Integration radius:%s px\n' "$RADIUS_VALUE"
printf '============================================================\n'

# ------------------------------------------------------------
# MAIN PIPELINE
# ------------------------------------------------------------

SENSOR_ENV="$SENSOR" \
CONSTANT_ENV="$CONSTANT" \
RADIUS_ENV="$RADIUS_VALUE" \
SOURCE_DATA_DIR_ENV="$SOURCE_DATA_DIR" \
TARGET_DATA_DIR_ENV="$TARGET_DATA_DIR" \
TARGET_MERGED_ROOT_ENV="$TARGET_MERGED_ROOT" \
SPOT_LUM_DIR_ENV="$SPOT_LUM_DIR" \
SPOT_LUM_MACRO_ENV="$SPOT_LUM_MACRO" \
COORD_MATCH_RADIUS_ENV="$COORD_MATCH_RADIUS" \
VBD_ENV="$VBD" \
"$PYTHON" <<'PY_PIPELINE'
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

sensor = os.environ["SENSOR_ENV"]
constant = os.environ["CONSTANT_ENV"]
radius = float(os.environ["RADIUS_ENV"])
source_data = Path(os.environ["SOURCE_DATA_DIR_ENV"])
target_data = Path(os.environ["TARGET_DATA_DIR_ENV"])
target_merged_root = Path(os.environ["TARGET_MERGED_ROOT_ENV"])
spot_lum_dir = Path(os.environ["SPOT_LUM_DIR_ENV"])
spot_lum_macro = os.environ["SPOT_LUM_MACRO_ENV"]
coord_match_radius = float(os.environ["COORD_MATCH_RADIUS_ENV"])
vbd_text = os.environ.get("VBD_ENV", "").strip()

integration_area = math.pi * radius * radius


# ============================================================
# HELPERS
# ============================================================

def phase_sort_key(phase: str):
    if phase == "before_annealing":
        return (0, -math.inf, -math.inf, phase)

    m = re.fullmatch(
        r"annealing_T=([-+]?(?:\d+(?:\.\d*)?|\.\d+))_h="
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+))",
        phase,
    )
    if m:
        return (1, float(m.group(1)), float(m.group(2)), phase)

    return (2, math.inf, math.inf, phase)


def extract_parameter_from_name(filename: str, parameter: str) -> float:
    m = re.search(
        rf"(?:^|_){re.escape(parameter)}="
        rf"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
        filename,
    )
    if m is None:
        raise ValueError(
            f"Cannot extract parameter {parameter!r} from filename: {filename}"
        )
    return float(m.group(1))


def ensure_link(link_path: Path, source_path: Path):
    """Create a non-destructive symlink to already-processed source data."""
    if not source_path.exists():
        return

    if link_path.is_symlink():
        try:
            if link_path.resolve() == source_path.resolve():
                return
        except FileNotFoundError:
            pass
        link_path.unlink()
    elif link_path.exists():
        # Never destroy a real user directory/file in the target tree.
        print(
            f"WARNING: {link_path} already exists and is not a symlink; "
            "leaving it untouched.",
            file=sys.stderr,
        )
        return

    link_path.symlink_to(source_path.resolve(), target_is_directory=source_path.is_dir())


def rewrite_coordinate_file(source: Path, target: Path):
    """
    Copy a 4-column coordinate file while imposing one arbitrary aperture.

    Expected fields:
        x, y, integration_area_pix2, integration_radius_pix

    x and y are copied exactly.  Both area and radius are replaced so that
    the metadata remain internally consistent with the requested aperture.
    """
    rows_written = 0
    target.parent.mkdir(parents=True, exist_ok=True)

    with source.open("r", encoding="utf-8") as src, target.open(
        "w", encoding="utf-8"
    ) as dst:
        for line_number, raw in enumerate(src, start=1):
            line = raw.strip()
            if not line:
                continue

            fields = [part.strip() for part in line.split(",")]
            if len(fields) < 4:
                raise ValueError(
                    f"Invalid coordinate line {line_number} in {source}: {raw.rstrip()}"
                )

            try:
                x = float(fields[0])
                y = float(fields[1])
                # Validate original numeric geometry even though it is replaced.
                float(fields[2])
                float(fields[3])
            except ValueError as exc:
                raise ValueError(
                    f"Non-numeric coordinate line {line_number} in {source}: "
                    f"{raw.rstrip()}"
                ) from exc

            dst.write(
                f"{x:.6f}, {y:.6f}, {integration_area:.6f}, {radius:.6f}\n"
            )
            rows_written += 1

    if rows_written == 0:
        raise RuntimeError(f"Coordinate file is empty: {source}")

    return rows_written


def nearest_one_to_one(result: pd.DataFrame, plan: pd.DataFrame, tolerance: float):
    """
    Associate ROOT output rows with the already-ID-labelled coordinate plan.

    IMPORTANT: this is NOT a global-ID assignment.  Global IDs already exist
    in 'plan'.  This only reconnects rows emitted by ROOT to the coordinates
    that were supplied to the macro in the same call.
    """
    pairs = []

    for r in result.itertuples(index=False):
        for p in plan.itertuples(index=False):
            d = math.hypot(float(r.x) - float(p.x), float(r.y) - float(p.y))
            if d <= tolerance:
                pairs.append((d, int(r.result_index), int(p.plan_index)))

    pairs.sort(key=lambda item: item[0])

    used_result = set()
    used_plan = set()
    mapping = {}

    for distance, result_index, plan_index in pairs:
        if result_index in used_result or plan_index in used_plan:
            continue
        mapping[result_index] = plan_index
        used_result.add(result_index)
        used_plan.add(plan_index)

    return mapping


# ============================================================
# DISCOVER SOURCE PHASE/RUN STRUCTURE
# ============================================================

phase_dirs = [
    p for p in source_data.iterdir()
    if p.is_dir()
    and (p.name == "before_annealing" or p.name.startswith("annealing_"))
]
phase_dirs.sort(key=lambda p: phase_sort_key(p.name))

run_pattern = re.compile(
    rf"^{re.escape(sensor)}_{re.escape(constant)}=(.+)_run=(.+)$"
)

datasets = []

for phase_dir in phase_dirs:
    matches = [
        p for p in phase_dir.iterdir()
        if p.is_dir() and run_pattern.fullmatch(p.name)
    ]

    if not matches:
        # The source analysis may legitimately not contain this sensor/condition
        # in every directory; only actual analysed phases are propagated.
        continue

    if len(matches) > 1:
        raise RuntimeError(
            f"More than one source run matches {sensor}_{constant}=*_run=* "
            f"inside {phase_dir}:\n  "
            + "\n  ".join(str(p) for p in matches)
        )

    run_dir = matches[0]
    m = run_pattern.fullmatch(run_dir.name)
    assert m is not None

    fixed_value = m.group(1)
    run_number = m.group(2)
    phase = phase_dir.name
    run_prefix = f"{sensor}_{constant}={fixed_value}"

    source_coord_dir = run_dir / "4coordinates"
    source_th2f_dir = run_dir / "5th2f"

    if not source_coord_dir.is_dir():
        raise RuntimeError(f"Missing source coordinate directory: {source_coord_dir}")
    if not source_th2f_dir.is_dir():
        raise RuntimeError(f"Missing source TH2F directory: {source_th2f_dir}")

    global_coord_files = sorted(
        source_coord_dir.glob("*_global_complete_coordinates.txt")
    )
    if len(global_coord_files) != 1:
        raise RuntimeError(
            f"Expected exactly one *_global_complete_coordinates.txt in "
            f"{source_coord_dir}; found {len(global_coord_files)}. "
            "Run the complete source pipeline through PHASE 3 first."
        )

    root_files = sorted(
        p for p in source_th2f_dir.glob("*data=diff*_processed_rotated_th2f.root")
        if "data=diffe" not in p.name
    )
    if not root_files:
        raise RuntimeError(f"No source data=diff TH2F files found in {source_th2f_dir}")

    datasets.append({
        "phase": phase,
        "phase_order": len(datasets),
        "source_phase_dir": phase_dir,
        "source_run_dir": run_dir,
        "fixed_value": fixed_value,
        "run_number": run_number,
        "run_prefix": run_prefix,
        "source_coord_dir": source_coord_dir,
        "source_global_coords": global_coord_files[0],
        "source_th2f_dir": source_th2f_dir,
        "root_files": root_files,
    })

if not datasets:
    raise RuntimeError(
        f"No source runs found for sensor={sensor}, constant={constant} in {source_data}."
    )

if datasets[0]["phase"] != "before_annealing":
    raise RuntimeError(
        "The source analysis must contain before_annealing for this sensor/condition."
    )

fixed_values = {d["fixed_value"] for d in datasets}
if len(fixed_values) != 1:
    details = "\n  ".join(
        f"{d['phase']}: {sensor}_{constant}={d['fixed_value']}"
        for d in datasets
    )
    raise RuntimeError(
        f"Inconsistent fixed {constant} values across source phases:\n  {details}"
    )

fixed_value = datasets[0]["fixed_value"]
condition_name = f"{sensor}_{constant}={fixed_value}"
source_merged_dir = source_data / "merged_files" / condition_name
target_merged_dir = target_merged_root / condition_name

if not source_merged_dir.is_dir():
    raise RuntimeError(
        "Source PHASE-3 merged directory not found:\n"
        f"  {source_merged_dir}\n"
        "Run the complete source pipeline first."
    )

# These files define the authoritative global-ID solution.  The new-radius
# pipeline must inherit them rather than solve the matching problem again.
source_geometry_plan_path = source_merged_dir / f"{sensor}_{constant}_geometry_plan.csv"
source_mapping_path = source_merged_dir / f"{sensor}_{constant}_global_mapping.csv"
source_catalog_path = source_merged_dir / f"{sensor}_{constant}_global_catalog.csv"
source_presence_path = source_merged_dir / f"{sensor}_{constant}_hotspot_presence.csv"
source_changes_path = source_merged_dir / f"{sensor}_{constant}_phase_changes.csv"

for required in (
    source_geometry_plan_path,
    source_mapping_path,
    source_catalog_path,
    source_presence_path,
    source_changes_path,
):
    if not required.is_file():
        raise RuntimeError(
            f"Required source PHASE-3 file not found: {required}"
        )

geometry_plan = pd.read_csv(source_geometry_plan_path)

required_plan_columns = {
    "spot", "phase", "run_number", "detection", "status", "x", "y"
}
missing = required_plan_columns.difference(geometry_plan.columns)
if missing:
    raise RuntimeError(
        f"Source geometry plan is missing columns {sorted(missing)}: "
        f"{source_geometry_plan_path}"
    )

# Reconstruct phase order from the source run structure; never infer it from
# lexical sorting of the CSV itself.
phase_to_order = {d["phase"]: d["phase_order"] for d in datasets}
unknown_phases = sorted(set(geometry_plan["phase"]) - set(phase_to_order))
if unknown_phases:
    raise RuntimeError(
        "Source geometry plan contains phases not found in the source run tree: "
        + ", ".join(unknown_phases)
    )

geometry_plan["phase_order"] = geometry_plan["phase"].map(phase_to_order).astype(int)
geometry_plan["spot"] = pd.to_numeric(geometry_plan["spot"], errors="raise").astype(int)
geometry_plan["x"] = pd.to_numeric(geometry_plan["x"], errors="raise")
geometry_plan["y"] = pd.to_numeric(geometry_plan["y"], errors="raise")

# Override the aperture globally while preserving EVERY global ID, coordinate,
# detection state, temporal status and matching diagnostic from the source.
geometry_plan["integration_radius"] = radius
geometry_plan["integration_area"] = integration_area
geometry_plan["integration_radius_source_phase"] = "user_fixed_radius"
geometry_plan["integration_radius_source_run"] = np.nan

# Fresh outputs only for this sensor/condition/radius. Other sensors already
# present in DATA_irradiated_isolated_R=<R> are left untouched.
if target_merged_dir.exists():
    shutil.rmtree(target_merged_dir)
target_merged_dir.mkdir(parents=True, exist_ok=True)

print(f"Source condition: {source_merged_dir}")
print(f"Target condition: {target_merged_dir}")
print(f"Inherited global hotspots: {geometry_plan['spot'].nunique()}")
print(f"Fixed integration radius: {radius:g} px")


# ============================================================
# CREATE TARGET TREE + RADIUS-SPECIFIC COORDINATE FILES
# ============================================================

for d in datasets:
    target_phase_dir = target_data / d["phase"]
    target_run_dir = target_phase_dir / d["source_run_dir"].name
    target_coord_dir = target_run_dir / "4coordinates"
    target_lum_dir = target_run_dir / "6luminosity"

    target_phase_dir.mkdir(parents=True, exist_ok=True)
    target_run_dir.mkdir(parents=True, exist_ok=True)

    # Keep the familiar run structure without copying large data or rerunning
    # preprocessing. These are read-only links to the already validated source.
    for dirname in ("1originals", "2processed", "3rotated", "5th2f"):
        ensure_link(target_run_dir / dirname, d["source_run_dir"] / dirname)

    # 4coordinates and 6luminosity are the radius-specific products and are
    # rebuilt on every execution for this selected run.
    if target_coord_dir.exists() or target_coord_dir.is_symlink():
        if target_coord_dir.is_symlink() or target_coord_dir.is_file():
            target_coord_dir.unlink()
        else:
            shutil.rmtree(target_coord_dir)
    if target_lum_dir.exists() or target_lum_dir.is_symlink():
        if target_lum_dir.is_symlink() or target_lum_dir.is_file():
            target_lum_dir.unlink()
        else:
            shutil.rmtree(target_lum_dir)

    target_coord_dir.mkdir(parents=True, exist_ok=True)
    target_lum_dir.mkdir(parents=True, exist_ok=True)

    # Copy every coordinate TXT so the target structure remains familiar.
    # All apertures are replaced by the requested fixed value.
    source_txt_files = sorted(d["source_coord_dir"].glob("*.txt"))
    if not source_txt_files:
        raise RuntimeError(f"No coordinate TXT files found in {d['source_coord_dir']}")

    copied_global = None
    for source_txt in source_txt_files:
        target_txt = target_coord_dir / source_txt.name
        rewrite_coordinate_file(source_txt, target_txt)
        if source_txt.name == d["source_global_coords"].name:
            copied_global = target_txt

    if copied_global is None or not copied_global.is_file():
        raise RuntimeError(
            f"Failed to create target global coordinate file for phase {d['phase']}."
        )

    d["target_run_dir"] = target_run_dir
    d["target_coord_dir"] = target_coord_dir
    d["target_lum_dir"] = target_lum_dir
    d["target_global_coords"] = copied_global

    print(
        f"Prepared {d['phase']} / {d['source_run_dir'].name}: "
        f"R={radius:g} px"
    )


# ============================================================
# ROOT LUMINOSITY RECALCULATION
# ============================================================

complete_phase_frames = []

for d in datasets:
    phase_plan = geometry_plan[
        geometry_plan["phase"] == d["phase"]
    ].copy()
    phase_plan = phase_plan.sort_values("spot", kind="stable").reset_index(drop=True)

    if phase_plan.empty:
        raise RuntimeError(f"No geometry-plan rows for phase {d['phase']}")

    # The source global_complete coordinate file is written sorted by global ID.
    # Verify the copied file has the same number of rows before calling ROOT.
    with d["target_global_coords"].open("r", encoding="utf-8") as handle:
        coord_count = sum(1 for line in handle if line.strip())
    if coord_count != len(phase_plan):
        raise RuntimeError(
            f"Coordinate/global-plan size mismatch in phase {d['phase']}: "
            f"{coord_count} coordinate rows vs {len(phase_plan)} global hotspots."
        )

    phase_measurements = []

    print(
        f"Recalculating {d['phase']}: {len(phase_plan)} global hotspots, "
        f"{len(d['root_files'])} TH2F files"
    )

    for root_file in d["root_files"]:
        t_value = extract_parameter_from_name(root_file.name, "T")
        v_value = extract_parameter_from_name(root_file.name, "v")

        result_csv = spot_lum_dir / "luminosity_results.csv"
        if result_csv.exists():
            result_csv.unlink()

        macro_call = (
            f'{spot_lum_macro}("{root_file}","{d["target_global_coords"]}")'
        )

        completed = subprocess.run(
            ["root", "-l", "-q", macro_call],
            cwd=spot_lum_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        if completed.returncode != 0:
            print(completed.stdout, file=sys.stderr)
            raise RuntimeError(
                f"ROOT luminosity macro failed for {root_file} "
                f"with exit code {completed.returncode}."
            )

        if not result_csv.is_file():
            print(completed.stdout, file=sys.stderr)
            raise RuntimeError(
                f"ROOT did not create luminosity_results.csv for {root_file}."
            )

        result = pd.read_csv(result_csv)
        result_csv.unlink()

        needed = {"x", "y", "luminosity", "error"}
        missing = needed.difference(result.columns)
        if missing:
            raise RuntimeError(
                f"ROOT output for {root_file} is missing columns {sorted(missing)}. "
                f"Found {list(result.columns)}"
            )

        for col in ("x", "y", "luminosity", "error"):
            result[col] = pd.to_numeric(result[col], errors="raise")

        result = result.reset_index(drop=True)
        result["result_index"] = np.arange(len(result), dtype=int)

        plan = phase_plan.reset_index(drop=True).copy()
        plan["plan_index"] = np.arange(len(plan), dtype=int)

        if len(result) != len(plan):
            raise RuntimeError(
                f"ROOT returned {len(result)} rows for {root_file}, but "
                f"{len(plan)} inherited global hotspots were supplied."
            )

        row_mapping = nearest_one_to_one(result, plan, coord_match_radius)
        if len(row_mapping) != len(result):
            raise RuntimeError(
                f"Could not reconnect every ROOT result to the inherited global-ID "
                f"coordinate plan for {root_file}. COORD_MATCH_RADIUS="
                f"{coord_match_radius} px."
            )

        plan_lookup = plan.set_index("plan_index")
        records = []

        for r in result.itertuples(index=False):
            p = plan_lookup.loc[row_mapping[int(r.result_index)]]

            record = {
                "spot": int(p["spot"]),
                "x": float(p["x"]),
                "y": float(p["y"]),
                "luminosity": float(r.luminosity),
                "error": float(r.error),
                "T": t_value,
                "v": v_value,
                "phase": p["phase"],
                "detection": p["detection"],
                "status": p["status"],
                "integration_radius": radius,
                "integration_area": integration_area,
                "run_number": d["run_number"],
            }

            # Preserve the original global-ID matching distance.  This is the
            # same diagnostic as in the authoritative source analysis.
            if "match_distance" in p.index:
                record["match_distance"] = p["match_distance"]

            # Preserve coordinate-propagation provenance when available.
            for col in (
                "geometry_source_phase",
                "geometry_source_run",
                "geometry_mode",
            ):
                if col in p.index:
                    record[col] = p[col]

            if constant == "T":
                record["v_fin"] = (
                    v_value - float(vbd_text) if vbd_text else np.nan
                )

            records.append(record)

        measured = pd.DataFrame(records)
        phase_measurements.append(measured)

        # Keep a per-operating-point CSV in the run's 6luminosity directory.
        stem = root_file.name.removesuffix("_th2f.root")
        point_output = d["target_lum_dir"] / f"luminosity_{stem}.csv"
        measured.sort_values("spot", kind="stable").to_csv(point_output, index=False)

    phase_complete = pd.concat(phase_measurements, ignore_index=True, sort=False)
    varying_column = "v" if constant == "T" else "T"
    phase_complete = phase_complete.sort_values(
        ["spot", varying_column], kind="stable"
    ).reset_index(drop=True)

    # Final per-run CSV in 6luminosity, matching the familiar source naming.
    run_output = d["target_lum_dir"] / (
        f"{d['run_prefix']}_{d['phase']}_run={d['run_number']}.csv"
    )
    phase_complete.to_csv(run_output, index=False)
    print(f"Created: {run_output}")

    # Phase-3 style per-phase output.
    phase_global_output = target_merged_dir / (
        f"{d['run_prefix']}_{d['phase']}_run={d['run_number']}_global_complete.csv"
    )
    phase_complete.to_csv(phase_global_output, index=False)
    print(f"Created: {phase_global_output}")

    phase_complete["__phase_order"] = d["phase_order"]
    complete_phase_frames.append(phase_complete)


# ============================================================
# PHASE 3 - MERGE WITH INHERITED GLOBAL IDS
# ============================================================

master = pd.concat(complete_phase_frames, ignore_index=True, sort=False)
varying_column = "v" if constant == "T" else "T"
master = master.sort_values(
    ["spot", "__phase_order", varying_column], kind="stable"
).reset_index(drop=True)
master = master.drop(columns=["__phase_order"])

master_output = target_merged_dir / f"{sensor}_{constant}_all_phases_global_ID.csv"
master.to_csv(master_output, index=False)
print(f"Created: {master_output}")

# Save a radius-specific geometry plan with the SAME global IDs/coordinates.
geometry_output = target_merged_dir / f"{sensor}_{constant}_geometry_plan.csv"
geometry_to_save = geometry_plan.drop(columns=["phase_order"]).copy()
geometry_to_save.to_csv(geometry_output, index=False)
print(f"Created: {geometry_output}")

# Copy the authoritative ID/detection diagnostics.  The global mapping's
# aperture metadata, if present, is overwritten with this requested radius so
# that every file inside the radius-specific folder is self-consistent.
mapping = pd.read_csv(source_mapping_path)
if "integration_radius" in mapping.columns:
    mapping["integration_radius"] = radius
if "integration_area" in mapping.columns:
    mapping["integration_area"] = integration_area
mapping_output = target_merged_dir / f"{sensor}_{constant}_global_mapping.csv"
mapping.to_csv(mapping_output, index=False)

for source_file, target_name in (
    (source_catalog_path, f"{sensor}_{constant}_global_catalog.csv"),
    (source_presence_path, f"{sensor}_{constant}_hotspot_presence.csv"),
    (source_changes_path, f"{sensor}_{constant}_phase_changes.csv"),
):
    shutil.copy2(source_file, target_merged_dir / target_name)


# ============================================================
# CROSS-CHECK: GLOBAL IDS MUST BE IDENTICAL TO SOURCE
# ============================================================

source_spots = set(pd.read_csv(source_catalog_path)["spot"].astype(int))
target_spots = set(master["spot"].astype(int))

if source_spots != target_spots:
    missing_in_target = sorted(source_spots - target_spots)
    unexpected = sorted(target_spots - source_spots)
    raise RuntimeError(
        "Global-ID invariant failed. This radius-specific run does not contain "
        "exactly the source global IDs.\n"
        f"Missing in target: {missing_in_target}\n"
        f"Unexpected in target: {unexpected}"
    )

# For every phase, the set of global IDs must also match the geometry plan.
for d in datasets:
    expected = set(
        geometry_plan.loc[geometry_plan["phase"] == d["phase"], "spot"].astype(int)
    )
    actual = set(
        master.loc[master["phase"] == d["phase"], "spot"].astype(int)
    )
    if expected != actual:
        raise RuntimeError(
            f"Global-ID phase invariant failed for {d['phase']}: "
            f"expected {len(expected)} IDs, found {len(actual)}."
        )

print("")
print("============================================")
print("FIXED-RADIUS RECALCULATION COMPLETED")
print(f"Sensor / condition: {condition_name}")
print(f"Integration radius: {radius:g} px")
print(f"Integration area:   {integration_area:.6f} px^2")
print(f"Phases processed:   {len(datasets)}")
print(f"Global hotspots:    {len(source_spots)}")
print("Global IDs:         IDENTICAL TO SOURCE (verified)")
print(f"Master file:        {master_output}")
print(f"Output directory:   {target_merged_dir}")
print("============================================")
PY_PIPELINE

printf '\n'
printf '============================================================\n'
printf 'PIPELINE COMPLETED SUCCESSFULLY\n'
printf 'Sensor:       %s\n' "$SENSOR"
printf 'Constant:     %s\n' "$CONSTANT"
printf 'Radius:       %s px\n' "$RADIUS_VALUE"
printf 'Output DATA:  %s\n' "$TARGET_DATA_DIR"
printf '============================================================\n'
