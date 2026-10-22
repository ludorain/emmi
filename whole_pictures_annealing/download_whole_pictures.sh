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
tail -n +2 "$CSV_FILE" | while IFS=';' read -r phase run sensor || [[ -n "$phase$run$sensor" ]]
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
    # COPY AND RENAME FILES INDEPENDENTLY
    # If one file is missing, keep the others.
    # --------------------------------------------------------

    remote_files=(
        "overlay.png"
        "stitched.data=light.tif"
        "stitched.data=diff-denoised.png"
        "stitched.data=diff-denoised.tif"
    )

    output_suffixes=(
        "overlay.png"
        "data=light.tif"
        "data=diff-denoised.png"
        "data=diff-denoised.tif"
    )

    downloaded=0
    missing=0

    for i in "${!remote_files[@]}"; do
        remote_file="${remote_files[$i]}"
        output_file="${output_dir}/${prefix}_run=${run}_${output_suffixes[$i]}"

        if sshpass -p "$PASSWORD" scp \
            "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_BASE}/${run}/${remote_file}" \
            "$tmp_dir/"; then

            if mv "$tmp_dir/$remote_file" "$output_file"; then
                echo "OK: $output_file"
                ((downloaded+=1))
            else
                echo "ERROR: failed to rename $remote_file for run $run"
                ((missing+=1))
            fi
        else
            echo "WARNING: unable to download $remote_file for run $run"
            ((missing+=1))
        fi
    done

    # --------------------------------------------------------
    # REMOVE TEMPORARY DIRECTORY
    # --------------------------------------------------------

    rm -rf "$tmp_dir"

    echo "Run $run: $downloaded downloaded, $missing missing or failed."

    echo "Done."
    echo

done


echo "============================================================"
echo "All runs processed."
echo "============================================================"