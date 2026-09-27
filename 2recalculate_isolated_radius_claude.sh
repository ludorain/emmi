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


def diagnose_coordinate_file(path: Path):
    """
    Re-parse a coordinate file EXACTLY the way spot_luminosity_sum_irradiated.C
    does: replace commas with spaces, then try to read four whitespace-
    separated numeric tokens per non-empty line. This mirrors the C++ macro's
    own parsing rule (a line that fails becomes 'Warning: riga non valida
    ignorata' inside ROOT and is silently dropped from its output), so any
    line flagged here is a line ROOT would also drop -- but here it is
    reported BEFORE calling ROOT, with the exact line number and content,
    instead of only showing up later as an unexplained row-count mismatch.
    """
    problems = []
    n_ok = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            tokens = line.replace(",", " ").split()
            if len(tokens) < 4:
                problems.append((line_number, raw.rstrip("\n"), "fewer than 4 tokens"))
                continue
            ok = True
            for tok in tokens[:4]:
                try:
                    float(tok)
                except ValueError:
                    ok = False
                    break
            if ok:
                n_ok += 1
            else:
                problems.append((line_number, raw.rstrip("\n"), "non-numeric token"))
    return n_ok, problems


def nearest_one_to_one(result: pd.DataFrame, plan: pd.DataFrame, tolerance: float):
    """
    Associate ROOT output rows with the already-ID-labelled coordinate plan.

    IMPORTANT: this is NOT a global-ID assignment.  Global IDs already exist
    in 'plan'.  This only reconnects rows emitted by ROOT to the coordinates
    that were supplied to the macro in the same call.

    NOTE: 'result' and 'plan' are NOT required to have the same length.  At a
    given fixed radius, ROOT may legitimately return fewer rows than the
    number of coordinates supplied (e.g. because a circle now overlaps a
    neighbour or exceeds the image edge at this radius, whereas it did not at
    the radius originally used to detect it). This function only produces a
    partial mapping in that case; the caller decides how to report the
    unmatched plan entries.
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


def reorder_complete_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Match the column layout of the main pipeline's global_complete files."""
    preferred = ["spot", "x", "y", "luminosity", "error", "T", "v"]
    if "v_fin" in df.columns:
        preferred.append("v_fin")
    preferred += [
        "phase", "detection", "status", "excluded_at_radius",
        "integration_radius", "integration_area",
        "geometry_source_phase", "geometry_source_constant",
        "geometry_source_run", "geometry_mode",
        "integration_radius_source_phase",
        "integration_radius_source_constant",
        "integration_radius_source_run",
        "run_number", "match_distance",
    ]
    if "cross_constant_only" in df.columns:
        preferred.append("cross_constant_only")

    missing = [c for c in preferred if c not in df.columns]
    if missing:
        raise RuntimeError(
            "Internal error: recalculated luminosity table is missing columns "
            f"required by the main-pipeline format: {missing}"
        )
    return df[preferred]


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

# These files define the authoritative solution produced by the main pipeline.
# The global catalog/mapping are SENSOR-LEVEL and are shared between T and v.
# The geometry plan, presence and phase-change reports remain condition-specific.
source_geometry_plan_path = source_merged_dir / f"{sensor}_{constant}_geometry_plan.csv"
source_mapping_path = source_merged_dir / f"{sensor}_{constant}_global_mapping.csv"
source_catalog_path = source_merged_dir / f"{sensor}_{constant}_global_catalog.csv"
source_presence_path = source_merged_dir / f"{sensor}_{constant}_hotspot_presence.csv"
source_changes_path = source_merged_dir / f"{sensor}_{constant}_phase_changes.csv"

source_shared_catalog_path = source_data / "merged_files" / f"{sensor}_global_catalog.csv"
source_shared_mapping_path = source_data / "merged_files" / f"{sensor}_global_mapping.csv"

for required in (
    source_geometry_plan_path,
    source_mapping_path,
    source_catalog_path,
    source_presence_path,
    source_changes_path,
    source_shared_catalog_path,
    source_shared_mapping_path,
):
    if not required.is_file():
        raise RuntimeError(
            f"Required source PHASE-3 file not found: {required}"
        )

geometry_plan = pd.read_csv(source_geometry_plan_path)

required_plan_columns = {
    "spot", "phase", "run_number", "constant",
    "detection", "status", "x", "y",
    "geometry_source_phase", "geometry_source_constant",
    "geometry_source_run", "geometry_mode",
    "integration_radius_source_phase",
    "integration_radius_source_constant",
    "integration_radius_source_run",
    "match_distance", "cross_constant_only",
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

plan_constants = set(geometry_plan["constant"].astype(str))
if plan_constants != {constant}:
    raise RuntimeError(
        f"Source geometry plan {source_geometry_plan_path} does not belong only "
        f"to constant={constant}. Found constants: {sorted(plan_constants)}"
    )

# Override ONLY the integration aperture.  Global IDs, coordinates, detection
# state/status and coordinate provenance remain exactly those of the main
# pipeline.  The aperture provenance is deliberately replaced, because this
# second pipeline imposes the radius by hand rather than deriving it from a
# detected hotspot.
geometry_plan["integration_radius"] = radius
geometry_plan["integration_area"] = integration_area
geometry_plan["integration_radius_source_phase"] = "user_fixed_radius"
geometry_plan["integration_radius_source_constant"] = "user_fixed_radius"
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

    # Preserve the source-finding coordinate files exactly as they were produced
    # by the main pipeline.  Their radii describe the original detected geometry
    # and must NOT be rewritten.  Only the global-complete coordinate file used
    # for the forced luminosity calculation receives the user-selected radius.
    source_txt_files = sorted(d["source_coord_dir"].glob("*.txt"))
    if not source_txt_files:
        raise RuntimeError(f"No coordinate TXT files found in {d['source_coord_dir']}")

    copied_global = None
    for source_txt in source_txt_files:
        target_txt = target_coord_dir / source_txt.name
        if source_txt.name == d["source_global_coords"].name:
            rewrite_coordinate_file(source_txt, target_txt)
            copied_global = target_txt
        else:
            shutil.copy2(source_txt, target_txt)

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

    # Catch, BEFORE calling ROOT, any line that ROOT's own parser would
    # silently reject (this is what previously surfaced only as an opaque
    # "ROOT returned N rows, but M were supplied" error after the fact).
    n_parseable, coord_problems = diagnose_coordinate_file(d["target_global_coords"])
    if coord_problems:
        print(
            f"WARNING: {len(coord_problems)}/{coord_count} line(s) in "
            f"{d['target_global_coords']} would NOT be read as valid "
            "coordinates by the ROOT macro's own parser (comma-split, then "
            "4 numeric tokens). These are almost certainly the rows that "
            "will end up missing from ROOT's output for every TH2F file in "
            f"this phase ({d['phase']}):",
            file=sys.stderr,
        )
        for line_number, content, reason in coord_problems[:20]:
            print(f"    line {line_number} ({reason}): {content!r}", file=sys.stderr)
        if len(coord_problems) > 20:
            print(f"    ... and {len(coord_problems) - 20} more.", file=sys.stderr)
    else:
        print(
            f"Coordinate file check OK: all {n_parseable} lines in "
            f"{d['target_global_coords'].name} parse as 4 valid numbers."
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

        # NOTE ON THE FIX (see chat explanation):
        # At a single fixed radius applied to EVERY inherited global hotspot,
        # ROOT can legitimately return FEWER rows than coordinates supplied
        # for a given file (typical causes: the enlarged circle now overlaps
        # a neighbouring hotspot, or exceeds the image edge, at this radius
        # even though it did not at the radius originally used to detect it).
        # This is expected behaviour of an "isolated-spot" analysis, not a
        # pipeline bug, so it must not abort the whole recalculation. Only a
        # ROOT output LARGER than what was supplied is treated as fatal,
        # since that can only mean a corrupted coordinate file or the wrong
        # macro/version being invoked.
        if len(result) > len(plan):
            raise RuntimeError(
                f"ROOT returned {len(result)} rows for {root_file}, MORE than "
                f"the {len(plan)} inherited global hotspots that were supplied. "
                "This cannot be reconciled automatically; check the coordinate "
                "file and the ROOT macro."
            )

        row_mapping = nearest_one_to_one(result, plan, coord_match_radius)
        if len(row_mapping) != len(result):
            raise RuntimeError(
                f"Could not reconnect every ROOT result to the inherited global-ID "
                f"coordinate plan for {root_file}. COORD_MATCH_RADIUS="
                f"{coord_match_radius} px."
            )

        matched_plan_indices = set(row_mapping.values())
        excluded_plan_indices = [
            idx for idx in plan["plan_index"] if idx not in matched_plan_indices
        ]

        if excluded_plan_indices:
            plan_lookup_for_warning = plan.set_index("plan_index")
            excluded_spots = sorted(
                int(plan_lookup_for_warning.loc[idx, "spot"])
                for idx in excluded_plan_indices
            )
            print(
                f"WARNING: {root_file.name}: ROOT returned {len(result)}/"
                f"{len(plan)} rows at radius={radius:g} px. Global spot(s) "
                "excluded for this file (likely overlapping a neighbour or "
                f"too close to the image edge at this radius): {excluded_spots}",
                file=sys.stderr,
            )

        plan_lookup = plan.set_index("plan_index")
        result_lookup = result.set_index("result_index")
        # Invert result_index -> plan_index into plan_index -> result_index
        # for a direct per-plan-row lookup below.
        plan_to_result = {
            plan_index: result_index
            for result_index, plan_index in row_mapping.items()
        }

        records = []

        # Iterate over every INHERITED plan row (not just the ROOT results),
        # so the output always keeps the full, stable set of global IDs for
        # this phase, matching the main pipeline's format. Spots ROOT could
        # not measure at this radius get NaN luminosity/error and are
        # flagged via 'excluded_at_radius', instead of silently shrinking
        # the table or aborting the run.
        for plan_index, p in plan_lookup.iterrows():
            result_index = plan_to_result.get(plan_index)
            matched = result_index is not None
            r = result_lookup.loc[result_index] if matched else None

            record = {
                "spot": int(p["spot"]),
                "x": float(p["x"]),
                "y": float(p["y"]),
                "luminosity": float(r["luminosity"]) if matched else float("nan"),
                "error": float(r["error"]) if matched else float("nan"),
                "T": t_value,
                "v": v_value,
                "phase": p["phase"],
                "detection": p["detection"],
                "status": p["status"],
                "excluded_at_radius": (not matched),
                "integration_radius": radius,
                "integration_area": integration_area,
                "run_number": d["run_number"],
            }

            # Preserve the original global-ID matching distance.  This is the
            # same diagnostic as in the authoritative source analysis.
            if "match_distance" in p.index:
                record["match_distance"] = p["match_distance"]

            # Preserve the complete coordinate-propagation provenance from the
            # main pipeline.  The integration-radius provenance in phase_plan has
            # already been replaced by user_fixed_radius above.
            for col in (
                "geometry_source_phase",
                "geometry_source_constant",
                "geometry_source_run",
                "geometry_mode",
                "integration_radius_source_phase",
                "integration_radius_source_constant",
                "integration_radius_source_run",
                "cross_constant_only",
            ):
                if col in p.index:
                    record[col] = p[col]

            if constant == "T":
                record["v_fin"] = (
                    v_value - float(vbd_text) if vbd_text else np.nan
                )

            records.append(record)

        measured = pd.DataFrame(records)
        measured = reorder_complete_columns(measured)
        phase_measurements.append(measured)

        # Keep a per-operating-point CSV in the run's 6luminosity directory,
        # using the same column layout as the main pipeline's complete outputs.
        stem = root_file.name.removesuffix("_th2f.root")
        point_output = d["target_lum_dir"] / f"luminosity_{stem}.csv"
        measured.sort_values("spot", kind="stable").to_csv(point_output, index=False)

    phase_complete = pd.concat(phase_measurements, ignore_index=True, sort=False)
    varying_column = "v" if constant == "T" else "T"
    phase_complete = phase_complete.sort_values(
        ["spot", varying_column], kind="stable"
    ).reset_index(drop=True)
    phase_complete = reorder_complete_columns(phase_complete)

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

    n_excluded_phase = int(phase_complete["excluded_at_radius"].sum())
    if n_excluded_phase:
        print(
            f"  -> {n_excluded_phase} (spot, operating point) row(s) in this "
            f"phase were excluded at radius={radius:g} px (see WARNINGs above)."
        )

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

# Copy the authoritative matching/detection diagnostics UNCHANGED.
# global_mapping describes how the source-finder detections were associated to
# sensor-level IDs; its original integration_radius/integration_area therefore
# belong to the detection geometry and must not be replaced by the arbitrary
# luminosity-integration radius.  The latter is recorded in geometry_plan and in
# every recalculated luminosity row.
shutil.copy2(
    source_mapping_path,
    target_merged_dir / f"{sensor}_{constant}_global_mapping.csv",
)

for source_file, target_name in (
    (source_catalog_path, f"{sensor}_{constant}_global_catalog.csv"),
    (source_presence_path, f"{sensor}_{constant}_hotspot_presence.csv"),
    (source_changes_path, f"{sensor}_{constant}_phase_changes.csv"),
):
    shutil.copy2(source_file, target_merged_dir / target_name)

# Also keep the SENSOR-LEVEL shared catalog and mapping in the radius-specific
# merged root.  These are the authoritative proof that T and v use the same
# global-ID namespace in the main pipeline.
shutil.copy2(
    source_shared_catalog_path,
    target_merged_root / f"{sensor}_global_catalog.csv",
)
shutil.copy2(
    source_shared_mapping_path,
    target_merged_root / f"{sensor}_global_mapping.csv",
)


# ============================================================
# CROSS-CHECK: GLOBAL IDS MUST BE IDENTICAL TO SOURCE
# ============================================================

# The SENSOR-LEVEL catalog is authoritative.  The per-condition catalog is kept
# only for backward compatibility by the main pipeline and must contain the same
# ID universe.
source_spots = set(pd.read_csv(source_shared_catalog_path)["spot"].astype(int))
condition_catalog_spots = set(pd.read_csv(source_catalog_path)["spot"].astype(int))
if condition_catalog_spots != source_spots:
    raise RuntimeError(
        "Source catalog inconsistency: the per-condition catalog does not match "
        "the sensor-level shared catalog. Rerun the main pipeline before this "
        "radius recalculation."
    )

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

n_excluded_total = int(master["excluded_at_radius"].sum())

print("")
print("============================================")
print("FIXED-RADIUS RECALCULATION COMPLETED")
print(f"Sensor / condition: {condition_name}")
print(f"Integration radius: {radius:g} px")
print(f"Integration area:   {integration_area:.6f} px^2")
print(f"Phases processed:   {len(datasets)}")
print(f"Global hotspots:    {len(source_spots)}")
print("Global IDs:         IDENTICAL TO SOURCE (verified)")
if n_excluded_total:
    print(
        f"NOTE: {n_excluded_total} (spot, operating point) row(s) across the "
        f"whole run have luminosity=NaN because ROOT excluded them at "
        f"radius={radius:g} px (see 'excluded_at_radius' column and the "
        "WARNING lines above)."
    )
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