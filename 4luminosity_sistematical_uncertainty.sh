#!/usr/bin/env bash

# ============================================================================
# 4luminosity_sistematical_uncertainty.sh
#
# Purpose
# -------
# Compute the symmetric systematic uncertainty on the luminosity measured
# with integration radius R=20, using R=16 and R=24 as variations:
#
#   deltaL1 = |L(R=24) - L(R=20)|
#   deltaL2 = |L(R=20) - L(R=16)|
#   deltaL  = max(deltaL1, deltaL2)
#
# Only the CSV files in DATA_irradiated_isolated_R=20 are modified.
# Before computing deltaL, the script verifies that R=16/20/24 rows refer to
# the same run and the same physical coordinates. This prevents stale radius
# datasets generated with an older global-ID convention from being mixed.
# All pre-existing columns are preserved exactly as strings, except the old
# deltaL_plus/deltaL_minus columns, which are removed if present. The columns
# deltaL1, deltaL2 and deltaL are then added/updated.
#
# Usage
# -----
#   ./4luminosity_sistematical_uncertainty.sh
#       -> all sensors (A1 A2 B1 B2), both constants (T and v)
#
#   ./4luminosity_sistematical_uncertainty.sh A1
#       -> sensor A1, both constants
#
#   ./4luminosity_sistematical_uncertainty.sh A1 T
#   ./4luminosity_sistematical_uncertainty.sh T A1
#       -> sensor A1, constant T
#
#   ./4luminosity_sistematical_uncertainty.sh T
#       -> all sensors, constant T
#
# Accepted sensor values: A1 A2 B1 B2
# Accepted constant values: T, v (case-insensitive)

# ============================================================================

set -u
set -o pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$SCRIPT_DIR"

R20_DIR="$BASE_DIR/DATA_irradiated_isolated_R=20"
R16_DIR="$BASE_DIR/DATA_irradiated_isolated_R=16"
R24_DIR="$BASE_DIR/DATA_irradiated_isolated_R=24"

T_VALUE="20"

ALL_SENSORS=(A1 A2 B1 B2)
ALL_CONSTANTS=(T v)

# ---------------------------------------------------------------------------
# Helper functions for command-line parsing
# ---------------------------------------------------------------------------
usage() {
    cat <<USAGE
Usage:
  $0
  $0 <sensor>
  $0 <constant>
  $0 <sensor> <constant>
  $0 <constant> <sensor>

Sensors:   A1 A2 B1 B2
Constants: T v T=20 v=5
USAGE
}

normalize_sensor() {
    local x
    x="$(printf '%s' "$1" | tr '[:lower:]' '[:upper:]')"
    case "$x" in
        A1|A2|B1|B2) printf '%s\n' "$x" ;;
        *) return 1 ;;
    esac
}

normalize_constant() {
    local x
    x="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')"
    case "$x" in
        t|t=20) printf '%s\n' "T" ;;
        v|v=5)  printf '%s\n' "v" ;;
        *) return 1 ;;
    esac
}

# ---------------------------------------------------------------------------
# Parse up to two positional arguments.
# The order is intentionally flexible to make the script robust to input use.
# ---------------------------------------------------------------------------
if (( $# > 2 )); then
    echo "ERROR: too many arguments." >&2
    usage >&2
    exit 2
fi

SELECTED_SENSOR=""
SELECTED_CONSTANT=""

for arg in "$@"; do
    if sensor="$(normalize_sensor "$arg" 2>/dev/null)"; then
        if [[ -n "$SELECTED_SENSOR" ]]; then
            echo "ERROR: more than one sensor was specified." >&2
            usage >&2
            exit 2
        fi
        SELECTED_SENSOR="$sensor"
        continue
    fi

    if constant="$(normalize_constant "$arg" 2>/dev/null)"; then
        if [[ -n "$SELECTED_CONSTANT" ]]; then
            echo "ERROR: more than one constant was specified." >&2
            usage >&2
            exit 2
        fi
        SELECTED_CONSTANT="$constant"
        continue
    fi

    echo "ERROR: unrecognized argument '$arg'." >&2
    usage >&2
    exit 2
done

if [[ -n "$SELECTED_SENSOR" ]]; then
    SENSORS=("$SELECTED_SENSOR")
else
    SENSORS=("${ALL_SENSORS[@]}")
fi

if [[ -n "$SELECTED_CONSTANT" ]]; then
    CONSTANTS=("$SELECTED_CONSTANT")
else
    CONSTANTS=("${ALL_CONSTANTS[@]}")
fi

# ---------------------------------------------------------------------------
# Verify that the three radius datasets exist before starting.
# ---------------------------------------------------------------------------
for d in "$R20_DIR" "$R16_DIR" "$R24_DIR"; do
    if [[ ! -d "$d" ]]; then
        echo "ERROR: required directory not found: $d" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# The delicate CSV work is done in embedded Python.
# Using the csv module (instead of pandas) avoids changing the textual values
# of pre-existing fields such as spot, luminosity, error, T and v.
# ---------------------------------------------------------------------------
python3 - "$BASE_DIR" "$T_VALUE" "${SENSORS[*]}" "${CONSTANTS[*]}" <<'PYCODE'
from __future__ import annotations

import csv
import os
import shutil
import sys
import tempfile
from decimal import Decimal, InvalidOperation
from pathlib import Path

BASE_DIR = Path(sys.argv[1])
T_VALUE = sys.argv[2]
SENSORS = sys.argv[3].split()
CONSTANTS = sys.argv[4].split()

RDIRS = {
    16: BASE_DIR / "DATA_irradiated_isolated_R=16",
    20: BASE_DIR / "DATA_irradiated_isolated_R=20",
    24: BASE_DIR / "DATA_irradiated_isolated_R=24",
}

NEW_COLUMNS = ["deltaL1", "deltaL2", "deltaL"]
DEPRECATED_COLUMNS = {"deltaL_plus", "deltaL_minus"}
REQUIRED_COLUMNS = {"spot", "x", "y", "luminosity", "T", "v", "phase", "run_number", "integration_radius"}

# Backups are created only once. Re-running the script therefore does not
# overwrite the original pre-systematics backup.
BACKUP_SUFFIX = ".pre_systematic_uncertainty.bak"

# Radius-specific files are generated from the same geometry plan, so their
# hotspot centres must agree to much better than one pixel.
COORD_TOL = Decimal("0.001")  # pixels
RADIUS_TOL = Decimal("0.000001")


class StructuralMismatch(RuntimeError):
    """The radius datasets do not share the same current ID/geometry solution."""


def validate_shared_catalogs(sensor: str) -> None:
    paths = {r: RDIRS[r] / "merged_files" / f"{sensor}_global_catalog.csv" for r in (16,20,24)}
    if not all(p.is_file() for p in paths.values()):
        missing = [str(p) for p in paths.values() if not p.is_file()]
        print(
            "WARNING: shared sensor-level catalog cross-check skipped because "
            f"these files are missing: {missing}"
        )
        return

    catalogs = {}
    for r, path in paths.items():
        with path.open("r", newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None or not {"spot","x_ref","y_ref"}.issubset(reader.fieldnames):
                raise StructuralMismatch(f"invalid shared catalog structure: {path}")
            cat = {}
            for row in reader:
                spot = normalize_numeric(row["spot"], "spot")
                cat[spot] = (Decimal(row["x_ref"].strip()), Decimal(row["y_ref"].strip()))
            catalogs[r] = cat

    ref = catalogs[20]
    for r in (16,24):
        if set(catalogs[r]) != set(ref):
            raise StructuralMismatch(
                f"sensor-level global-ID catalog differs between R=20 and R={r}. "
                "Regenerate all radius datasets from the same current primary pipeline."
            )
        for spot, (x20,y20) in ref.items():
            xr, yr = catalogs[r][spot]
            if max(abs(xr-x20), abs(yr-y20)) > COORD_TOL:
                raise StructuralMismatch(
                    f"sensor-level catalog geometry mismatch for sensor={sensor}, spot={spot}: "
                    f"R20=({x20},{y20}), R{r}=({xr},{yr}). Regenerate R16/R20/R24."
                )


def condition_folder(sensor: str, constant: str) -> str:
    if constant == "T":
        return f"{sensor}_T={T_VALUE}"

    if constant == "v":
        V_VALUES = {
            "A1": "5",
            "A2": "5",
            "B1": "3",
            "B2": "3",
        }

        if sensor not in V_VALUES:
            raise ValueError(f"Unsupported sensor: {sensor}")

        return f"{sensor}_v={V_VALUES[sensor]}"

    raise ValueError(f"Unsupported constant: {constant}")


def normalize_numeric(text: str, field: str) -> str:
    """Normalize numeric key fields without changing what is written to CSV."""
    raw = text.strip()
    try:
        value = Decimal(raw)
    except InvalidOperation as exc:
        raise ValueError(f"Field '{field}' is not numeric: {text!r}") from exc

    # Decimal('20.0') and Decimal('20') must match.
    if value == 0:
        return "0"
    return str(value.normalize())


def row_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    """Identity of one measurement across the three radius datasets."""
    return (
        normalize_numeric(row["spot"], "spot"),
        row["phase"].strip(),
        normalize_numeric(row["T"], "T"),
        normalize_numeric(row["v"], "v"),
        row.get("run_number", "").strip(),
    )


def all_phases_key(row: dict[str, str]) -> tuple[str, ...]:
    """Use exactly the same measurement identity in run and master files."""
    return row_key(row)


def decimal_to_text(value: Decimal) -> str:
    """Stable text representation for newly computed uncertainty columns."""
    if value == 0:
        return "0"
    # Fixed-point avoids unnecessary scientific notation for ordinary values.
    text = format(value.normalize(), "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("CSV has no header")
        fieldnames = list(reader.fieldnames)
        missing = REQUIRED_COLUMNS.difference(fieldnames)
        if missing:
            raise ValueError(f"missing required columns: {sorted(missing)}")
        rows = list(reader)
    return fieldnames, rows


def build_unique_index(rows: list[dict[str, str]], path: Path) -> dict[tuple[str, ...], dict[str, str]]:
    index: dict[tuple[str, ...], dict[str, str]] = {}
    duplicates: list[tuple[str, ...]] = []

    for row in rows:
        key = row_key(row)
        if key in index:
            duplicates.append(key)
        else:
            index[key] = row

    if duplicates:
        preview = ", ".join(map(str, duplicates[:5]))
        raise ValueError(
            f"non-unique matching key (spot, phase, T, v, run_number) in {path.name}. "
            f"Examples: {preview}"
        )
    return index


def atomic_write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    """Write to a temporary file and atomically replace the original."""
    backup = Path(str(path) + BACKUP_SUFFIX)
    if not backup.exists():
        shutil.copy2(path, backup)

    fd, temp_name = tempfile.mkstemp(
        prefix=path.name + ".tmp.", dir=str(path.parent), text=True
    )
    os.close(fd)
    temp_path = Path(temp_name)

    try:
        with temp_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="raise")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temp_path, path)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def as_decimal(row: dict[str, str], field: str) -> Decimal:
    try:
        return Decimal(row[field].strip())
    except (InvalidOperation, KeyError) as exc:
        raise ValueError(f"invalid numeric field {field!r}: {row.get(field)!r}") from exc


def validate_same_measurement(row16: dict[str, str], row20: dict[str, str],
                              row24: dict[str, str], key: tuple[str, ...]) -> None:
    """Reject stale/mixed radius files before luminosities are compared."""
    # The centres must be the same because all fixed-radius recalculations are
    # required to inherit one common geometry plan from the primary pipeline.
    x16, x20, x24 = (as_decimal(r, "x") for r in (row16, row20, row24))
    y16, y20, y24 = (as_decimal(r, "y") for r in (row16, row20, row24))
    if max(abs(x16-x20), abs(x24-x20), abs(y16-y20), abs(y24-y20)) > COORD_TOL:
        raise StructuralMismatch(
            f"geometry mismatch for key {key}: "
            f"R16=({x16},{y16}), R20=({x20},{y20}), R24=({x24},{y24}). "
            "The radius datasets were probably produced with different global-ID/geometry solutions; "
            "regenerate R16, R20 and R24 from the same current source pipeline."
        )

    expected = ((row16, Decimal("16")), (row20, Decimal("20")), (row24, Decimal("24")))
    for row, target in expected:
        actual = as_decimal(row, "integration_radius")
        if abs(actual-target) > RADIUS_TOL:
            raise StructuralMismatch(
                f"wrong integration_radius for key {key}: expected R={target}, found {actual}"
            )

    # Provenance/detection metadata should also be identical across radius copies.
    for field in (
        "detection", "status", "geometry_source_phase", "geometry_source_constant",
        "geometry_source_run", "geometry_mode", "cross_constant_only"
    ):
        vals = [r.get(field, "").strip() for r in (row16, row20, row24)]
        nonempty = [v for v in vals if v != ""]
        if nonempty and len(set(nonempty)) != 1:
            raise StructuralMismatch(f"metadata mismatch in field {field!r} for key {key}: {vals}")


def process_one_run_file(path20: Path, path16: Path, path24: Path):
    """
    Compute systematic uncertainties for a single *_global_complete.csv file.

    Safety policy:
      - if either reference file is missing -> skip this file;
      - if keys are duplicated -> skip this file;
      - if a R=20 row has no exact R=16 or R=24 match -> skip this file;
      - only after all checks pass is the R=20 file modified.
    """
    if not path16.is_file() or not path24.is_file():
        missing = []
        if not path16.is_file():
            missing.append("R=16")
        if not path24.is_file():
            missing.append("R=24")
        raise FileNotFoundError("missing corresponding file in " + " and ".join(missing))

    fields20, rows20 = read_csv(path20)
    _, rows16 = read_csv(path16)
    _, rows24 = read_csv(path24)

    idx16 = build_unique_index(rows16, path16)
    idx20 = build_unique_index(rows20, path20)
    idx24 = build_unique_index(rows24, path24)

    # Do not require R16/R24 to have exactly the same full key set: extra rows
    # are harmless. Every R20 row, however, must have both matching references.
    missing16 = [k for k in idx20 if k not in idx16]
    missing24 = [k for k in idx20 if k not in idx24]
    if missing16 or missing24:
        pieces = []
        if missing16:
            pieces.append(f"{len(missing16)} R=20 keys missing in R=16")
        if missing24:
            pieces.append(f"{len(missing24)} R=20 keys missing in R=24")
        raise ValueError("; ".join(pieces))

    uncertainty_by_all_key: dict[tuple[str, ...], tuple[str, str, str]] = {}
    diagnostics: list[dict[str, str]] = []

    for row20 in rows20:
        key = row_key(row20)
        row16 = idx16[key]
        row24 = idx24[key]

        validate_same_measurement(row16, row20, row24, key)

        try:
            L16 = Decimal(row16["luminosity"].strip())
            L20 = Decimal(row20["luminosity"].strip())
            L24 = Decimal(row24["luminosity"].strip())
        except InvalidOperation as exc:
            raise ValueError(f"invalid luminosity for key {key}") from exc

        delta1 = abs(L24 - L20)
        delta2 = abs(L20 - L16)
        delta = max(delta1, delta2)

        diagnostics.append({
            "file": path20.name,
            "spot": row20["spot"],
            "phase": row20["phase"],
            "run_number": row20.get("run_number", ""),
            "T": row20["T"],
            "v": row20["v"],
            "x": row20["x"],
            "y": row20["y"],
            "detection": row20.get("detection", ""),
            "L16": decimal_to_text(L16),
            "L20": decimal_to_text(L20),
            "L24": decimal_to_text(L24),
            "deltaL1": decimal_to_text(delta1),
            "deltaL2": decimal_to_text(delta2),
            "deltaL": decimal_to_text(delta),
            "deltaL_over_absL20": (
                decimal_to_text(delta / abs(L20)) if L20 != 0 else ""
            ),
        })

        row20["deltaL1"] = decimal_to_text(delta1)
        row20["deltaL2"] = decimal_to_text(delta2)
        row20["deltaL"] = decimal_to_text(delta)

        # Remove deprecated asymmetric columns from the row if this script
        # is re-run on CSV files produced by an older pipeline version.
        for old_col in DEPRECATED_COLUMNS:
            row20.pop(old_col, None)

        # The same full key, including run_number, is propagated to the all-phases file.
        akey = all_phases_key(row20)
        values = (row20["deltaL1"], row20["deltaL2"], row20["deltaL"])
        old = uncertainty_by_all_key.get(akey)
        if old is not None and old != values:
            raise ValueError(
                f"conflicting systematic uncertainties for all-phases key {akey}"
            )
        uncertainty_by_all_key[akey] = values

    # Remove the deprecated asymmetric columns from the output header.
    out_fields = [col for col in fields20 if col not in DEPRECATED_COLUMNS]
    for col in NEW_COLUMNS:
        if col not in out_fields:
            out_fields.append(col)

    atomic_write_csv(path20, out_fields, rows20)
    return len(rows20), uncertainty_by_all_key, diagnostics


def update_all_phases_file(path: Path, mapping: dict[tuple[str, ...], tuple[str, str, str]]) -> tuple[int, int]:
    """Update an existing *_all_phases_global_ID.csv without recreating it."""
    fields, rows = read_csv(path)

    # Check uniqueness before modifying anything.
    seen: set[tuple[str, ...]] = set()
    duplicates: list[tuple[str, ...]] = []
    for row in rows:
        key = all_phases_key(row)
        if key in seen:
            duplicates.append(key)
        seen.add(key)
    if duplicates:
        preview = ", ".join(map(str, duplicates[:5]))
        raise ValueError(
            "non-unique key in all-phases file even after run_number "
            f"disambiguation. Examples: {preview}"
        )

    matched = 0
    unmatched = 0
    for row in rows:
        key = all_phases_key(row)
        values = mapping.get(key)
        if values is None:
            # Leave the new systematic values untouched if the row was not part
            # of the currently processed run files. This allows processing
            # subsets of sensors/constants without destroying earlier results.
            unmatched += 1
        else:
            row["deltaL1"], row["deltaL2"], row["deltaL"] = values
            matched += 1

        # Always remove deprecated asymmetric columns from the all-phases file.
        for old_col in DEPRECATED_COLUMNS:
            row.pop(old_col, None)

    out_fields = [col for col in fields if col not in DEPRECATED_COLUMNS]
    for col in NEW_COLUMNS:
        if col not in out_fields:
            out_fields.append(col)

    atomic_write_csv(path, out_fields, rows)
    return matched, unmatched


def find_all_phases_files(folder: Path, sensor: str, constant: str) -> list[Path]:
    # The expected names are e.g. A1_T_all_phases_global_ID.csv and
    # A1_v_all_phases_global_ID.csv. The slightly broader pattern tolerates
    # harmless naming additions while staying inside the selected folder.
    prefix = f"{sensor}_{constant}"
    return sorted(folder.glob(f"{prefix}*all_phases_global_ID.csv"))


total_run_files = 0
total_rows = 0
total_skipped = 0
total_all_phase_files = 0

print("=== Systematic luminosity uncertainty ===")
print(f"Base directory: {BASE_DIR}")
print(f"Sensors: {', '.join(SENSORS)}")
print(f"Constants: {', '.join(CONSTANTS)}")
print()

for sensor in SENSORS:
    validate_shared_catalogs(sensor)
    print(f"Shared global-ID catalog check passed for {sensor}.")
    for constant in CONSTANTS:
        folder_name = condition_folder(sensor, constant)
        folder20 = RDIRS[20] / "merged_files" / folder_name
        folder16 = RDIRS[16] / "merged_files" / folder_name
        folder24 = RDIRS[24] / "merged_files" / folder_name

        print(f"--- {sensor}, constant {constant} -> {folder_name} ---")

        if not folder20.is_dir():
            print(f"WARNING: R=20 folder not found, skipping: {folder20}")
            print()
            continue

        # Only phase/run files are processed here. The all-phases file is
        # updated separately from the resulting uncertainty map.
        run_files = sorted(folder20.glob("*_global_complete.csv"))
        if not run_files:
            print("WARNING: no *_global_complete.csv files found; nothing to do.")
            print()
            continue

        global_mapping: dict[tuple[str, ...], tuple[str, str, str]] = {}
        successful_files = 0
        condition_diagnostics: list[dict[str, str]] = []

        for path20 in run_files:
            path16 = folder16 / path20.name
            path24 = folder24 / path20.name

            try:
                nrows, local_mapping, local_diagnostics = process_one_run_file(path20, path16, path24)
            except StructuralMismatch as exc:
                print(f"FATAL structural mismatch: {path20.name}", file=sys.stderr)
                print(f"  {exc}", file=sys.stderr)
                raise
            except Exception as exc:
                total_skipped += 1
                print(f"SKIP: {path20.name}")
                print(f"      reason: {exc}")
                continue

            # Merge mappings conservatively. Conflicts indicate a structural
            # ambiguity and are never silently overwritten.
            conflict = False
            for key, pair in local_mapping.items():
                old = global_mapping.get(key)
                if old is not None and old != pair:
                    print(f"ERROR: conflicting all-phases mapping for key {key}")
                    conflict = True
                    break
            if conflict:
                # The run file itself has already been safely updated, but the
                # ambiguous row is not propagated to the all-phases file.
                continue

            global_mapping.update(local_mapping)
            condition_diagnostics.extend(local_diagnostics)
            successful_files += 1
            total_run_files += 1
            total_rows += nrows
            print(f"OK:   {path20.name} ({nrows} rows)")

        # Update the already-existing all-phases CSV rather than recreating it.
        # This minimizes changes to ordering and to all pre-existing values.
        if global_mapping:
            all_phase_files = find_all_phases_files(folder20, sensor, constant)
            if not all_phase_files:
                print("WARNING: no *_all_phases_global_ID.csv file found.")
            else:
                for all_path in all_phase_files:
                    try:
                        matched, unmatched = update_all_phases_file(all_path, global_mapping)
                    except Exception as exc:
                        print(f"SKIP all-phases: {all_path.name}")
                        print(f"                 reason: {exc}")
                        continue
                    total_all_phase_files += 1
                    print(
                        f"OK all-phases: {all_path.name} "
                        f"({matched} rows updated, {unmatched} rows unchanged)"
                    )

        if condition_diagnostics:
            diagnostic_path = folder20 / f"{sensor}_{constant}_systematic_uncertainty_diagnostics.csv"
            diagnostic_fields = [
                "file", "spot", "phase", "run_number", "T", "v", "x", "y",
                "detection", "L16", "L20", "L24", "deltaL1", "deltaL2", "deltaL",
                "deltaL_over_absL20"
            ]
            with diagnostic_path.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=diagnostic_fields)
                writer.writeheader()
                writer.writerows(condition_diagnostics)
            print(f"Diagnostic table: {diagnostic_path}")

        print(f"Completed {successful_files}/{len(run_files)} run files in {folder_name}.")
        print()

print("=== Final summary ===")
print(f"Run files updated:       {total_run_files}")
print(f"Rows updated:            {total_rows}")
print(f"Run files skipped:       {total_skipped}")
print(f"All-phases files updated:{total_all_phase_files}")
print(f"Backups suffix:          {BACKUP_SUFFIX}")
PYCODE

status=$?
if (( status != 0 )); then
    echo "ERROR: systematic-uncertainty procedure terminated with status $status." >&2
    exit "$status"
fi

echo "Done."
