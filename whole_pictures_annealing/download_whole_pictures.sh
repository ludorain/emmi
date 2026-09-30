#!/usr/bin/env bash

set -u
set -o pipefail

# ============================================================
# Download whole-device EMMI pictures for multiple runs
#
# Execute this script from:
#   whole_pictures_annealing/
#
# Expected structure:
#
#   whole_pictures_annealing/
#   ├── download_whole_pictures.sh
#   ├── runs.csv
#   ├── A1/
#   ├── A2/
#   ├── B1/
#   └── B2/
#
# CSV format:
#
# Phase,Run number,Sensor
# bef_ann,20260512-123812,A1
# annealing_T=75_h=5,20260513-123456,A1
# ...
#
# ============================================================


# ============================================================
# USER SETTINGS
# ============================================================

CSV_FILE="download_files_whole.csv"

REMOTE_USER="eic"
REMOTE_HOST="pccerna"
REMOTE_BASE="/home/eic/DATA/EMMI/actual"

PASSWORD="sipm4all"


# ============================================================
# CHECKS
# ============================================================

if ! command -v sshpass >/dev/null 2>&1; then
    echo "ERROR: sshpass is not installed."
    echo
    echo "On macOS with Homebrew you can install it with:"
    echo "  brew install hudochenkov/sshpass/sshpass"
    exit 1
fi

if [[ ! -f "$CSV_FILE" ]]; then
    echo "ERROR: CSV file not found: $CSV_FILE"
    exit 1
fi


# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

mkdir -p A1 A2 B1 B2


# ============================================================
# READ CSV
# ============================================================

echo
echo "============================================================"
echo "Starting download"
echo "CSV: $CSV_FILE"
echo "============================================================"
echo


# Skip first line because it contains the header.
tail -n +2 "$CSV_FILE" | while IFS=';' read -r phase run sensor
do

    # Remove possible Windows carriage returns and surrounding spaces.
    phase=$(echo "$phase" | sed 's/\r//g' | xargs)
    run=$(echo "$run" | sed 's/\r//g' | xargs)
    sensor=$(echo "$sensor" | sed 's/\r//g' | xargs)


    # --------------------------------------------------------
    # Skip empty lines
    # --------------------------------------------------------

    if [[ -z "$phase" && -z "$run" && -z "$sensor" ]]; then
        continue
    fi


    # --------------------------------------------------------
    # Validate sensor
    # --------------------------------------------------------

    case "$sensor" in
        A1|A2|B1|B2)
            ;;
        *)
            echo "WARNING: unknown sensor '$sensor'. Skipping run $run."
            echo
            continue
            ;;
    esac


    # --------------------------------------------------------
    # Build prefix
    # --------------------------------------------------------

    prefix="${sensor}_${phase}_whole"

    output_dir="$sensor"


    echo "------------------------------------------------------------"
    echo "Sensor : $sensor"
    echo "Phase  : $phase"
    echo "Run    : $run"
    echo "Prefix : $prefix"
    echo "Output : $output_dir/"
    echo "------------------------------------------------------------"


    # --------------------------------------------------------
    # Temporary directory
    #
    # This avoids problems if files with the same original
    # names are downloaded for multiple runs.
    # --------------------------------------------------------

    tmp_dir=$(mktemp -d)


    # --------------------------------------------------------
    # COPY FILES
    # --------------------------------------------------------

    sshpass -p "$PASSWORD" scp \
        "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_BASE}/${run}/overlay.png" \
        "$tmp_dir/"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: failed to download overlay.png for run $run"
        rm -rf "$tmp_dir"
        continue
    fi

    sshpass -p "$PASSWORD" scp \
        "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_BASE}/${run}/stitched.data=light.tif" \
        "$tmp_dir/"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: failed to download stitched.data=light.tif for run $run"
        rm -rf "$tmp_dir"
        continue
    fi


    # Check whether scp succeeded.
    if [[ $? -ne 0 ]]; then
        echo "ERROR: scp failed for run $run."
        rm -rf "$tmp_dir"
        echo
        continue
    fi


    # --------------------------------------------------------
    # RENAME FILES
    # --------------------------------------------------------

    mv \
        "$tmp_dir/overlay.png" \
        "${output_dir}/${prefix}_run=${run}_overlay.png"

    mv \
        "$tmp_dir/stitched.data=light.tif" \
        "${output_dir}/${prefix}_run=${run}_data=light.tif"


    # --------------------------------------------------------
    # REMOVE TEMPORARY DIRECTORY
    # --------------------------------------------------------

    rm -rf "$tmp_dir"


    echo "Done."
    echo

done


echo "============================================================"
echo "All runs processed."
echo "============================================================"