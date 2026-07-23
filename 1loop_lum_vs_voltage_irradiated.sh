#!/bin/bash
#ricorda che per A1: v_fin = v - 51.3, per B1: v_fin = v - 50.9

# --- CONFIGURAZIONE PERCORSI ---
BASE_DIR=$(pwd)
# Modificato per leggere dalla cartella T=20
ORIGINALS_DIR="$BASE_DIR/DATA_irradiated/annealing_T=100_h=25/A1_T=20_20260625-095817/1originals"
PROCESSED_DIR="$BASE_DIR/DATA_irradiated/annealing_T=100_h=25/A1_T=20_20260625-095817/2processed"
COORD_DIR="$BASE_DIR/DATA_irradiated/annealing_T=100_h=25/A1_T=20_20260625-095817/3coordinates"
ROOT_DIR="$BASE_DIR/DATA_irradiated/annealing_T=100_h=25/A1_T=20_20260625-095817/4th2f"
LUM_DIR="$BASE_DIR/DATA_irradiated/annealing_T=100_h=25/A1_T=20_20260625-095817/luminosity"

mkdir -p "$PROCESSED_DIR" "$COORD_DIR" "$ROOT_DIR" "$LUM_DIR"

# --- FASE 1: TROVARE IL FILE CON IL SECONDO VALORE MASSIMO DI V ---
echo "Searching for the file with the second-highest v..."

reference_file=""
reference_coords="$COORD_DIR/reference_coordinates.txt"

# Estrazione prefisso comune
first_file=$(find "$ORIGINALS_DIR" -maxdepth 1 -name '*_data=diff.tif' | head -n 1)

if [ -z "$first_file" ]; then
    echo "Error: no *_data=diff.tif files found in:"
    echo "$ORIGINALS_DIR"
    exit 1
fi

common_prefix=$(basename "$first_file" | cut -d'_' -f1)

# Crea una lista "v <TAB> percorso", la ordina per v decrescente,
# elimina eventuali valori di v duplicati e seleziona il secondo.
second_entry=$(
    for f in "$ORIGINALS_DIR"/*_data=diff.tif; do
        filename=$(basename "$f")
        val_v=$(echo "$filename" | sed -n 's/.*v=\([0-9.]*\).*/\1/p')

        if [ -n "$val_v" ]; then
            printf '%s\t%s\n' "$val_v" "$f"
        fi
    done |
    sort -t $'\t' -k1,1gr |
    awk -F '\t' '!seen[$1]++' |
    sed -n '2p'
)

if [ -z "$second_entry" ]; then
    echo "Error: fewer than two distinct v values were found."
    exit 1
fi

IFS=$'\t' read -r second_max_v reference_file <<< "$second_entry"

echo "Second-highest v found: $second_max_v"
echo "Reference file: $(basename "$reference_file")"
echo "Prefix: $common_prefix"


# --- FASE 2: GENERAZIONE COORDINATE DI RIFERIMENTO ---
ref_basename=$(basename "$reference_file" .tif)
ref_processed="$PROCESSED_DIR/${ref_basename}_processed.tif"

if [ ! -f "$ref_processed" ]; then
    cd "$BASE_DIR/manipulate_images/" || exit 1

    python process_image.py \
        --input "$reference_file" \
        --process remove_column_bias remove_hot_pixels remove_cold_pixels \
        --output "$ref_processed"
fi

cd "$BASE_DIR/find_centers/" || exit 1

python find_centers_irradiated.py \
    --input "$ref_processed" \
    --convolution \
    --coordinates_root "$reference_coords"

# --- FASE 3: LOOP DI PROCESSAMENTO SU TUTTI I FILE ---
cd "$BASE_DIR"
for input_path in "$ORIGINALS_DIR"/*_data=diff.tif; do
    filename=$(basename "$input_path")
    basename="${filename%.tif}"
    error_file="${input_path%diff.tif}diffe.tif"
    
    if [ ! -f "$error_file" ]; then continue; fi

    val_T=$(echo "$filename" | sed -n 's/.*T=\([0-9.]*\).*/\1/p')
    val_v=$(echo "$filename" | sed -n 's/.*v=\([0-9.]*\).*/\1/p')

    # Processamento
    if [ ! -f "$PROCESSED_DIR/${basename}_processed.tif" ]; then
        cd "$BASE_DIR/manipulate_images/"
        python process_image.py --input "$input_path" --process remove_column_bias remove_hot_pixels remove_cold_pixels --output "$PROCESSED_DIR/${basename}_processed.tif"
    fi
    
    cd "$BASE_DIR/manipulate_images/"
    python tif2th2.py --input "$PROCESSED_DIR/${basename}_processed.tif" --error "$error_file" --output "$ROOT_DIR/${basename}_processed_th2f.root"

    cd "$BASE_DIR/spot_luminosity/"
    root -l -q "spot_luminosity_sum_irradiated.C(\"$ROOT_DIR/${basename}_processed_th2f.root\",\"$reference_coords\")"

    if [ -f "luminosity_results.csv" ]; then
        [[ -z "$val_T" || -z "$val_v" ]] && { val_T="0"; val_v="0"; }
        # Manteniamo T e v originali nel CSV temporaneo
        awk -v t="$val_T" -v v="$val_v" 'BEGIN { FS=","; OFS="," } { gsub(/[[:space:]]+$/, "", $0); if ($0 != "") { if (NR == 1) print $0, "T", "v"; else print $0, t, v; } }' luminosity_results.csv > "$LUM_DIR/luminosity_${basename}_processed.csv"
        rm -f luminosity_results.csv
    fi
done

# --- FASE 4: UNIONE, CALCOLO V_FIN E ORDINAMENTO ---
echo "Merging and calculating v_fin..."
cd "$LUM_DIR"

python3 <<EOF
import pandas as pd
import glob
import os

# 1. Unione CSV
files = glob.glob("luminosity_*_processed.csv")
df_list = []
for f in files:
    df_list.append(pd.read_csv(f))

full_df = pd.concat(df_list, ignore_index=True)

# 2. Assegnazione ID spot
full_df['spot'] = full_df.groupby(['x', 'y']).ngroup()
cols = ['spot'] + [c for c in full_df.columns if c != 'spot']
full_df = full_df[cols]

# 3. Calcolo v_fin = v - 51.3
# Convertiamo in numerico per sicurezza
full_df['v'] = pd.to_numeric(full_df['v'], errors='coerce')
full_df['T'] = pd.to_numeric(full_df['T'], errors='coerce')
full_df['v_fin'] = full_df['v'] - 51.3

# 4. Ordinamento per spot e per v
full_df = full_df.sort_values(by=['spot', 'v'])

# 5. Salvataggio finale
output_name = "${common_prefix}_total.csv"
full_df.to_csv(output_name, index=False)
print(f"Final file created: {output_name}")
EOF

echo "Process completed successfully."