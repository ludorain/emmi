#!/bin/bash

BASE_DIR=$(pwd)
ORIGINALS_DIR="$BASE_DIR/files_v=9/1originals"
PROCESSED_DIR="$BASE_DIR/files_v=9/2processed"
COORD_DIR="$BASE_DIR/files_v=9/3coordinates"
ROOT_DIR="$BASE_DIR/files_v=9/4th2f"
LUM_DIR="$BASE_DIR/files_v=9/luminosity"

mkdir -p "$PROCESSED_DIR" "$COORD_DIR" "$ROOT_DIR" "$LUM_DIR"

# --- FASE 1: TROVARE IL FILE CON T MASSIMO ---
echo "Searching for the file with maximum T..."
max_T=-1
reference_file=""
reference_coords="$COORD_DIR/reference_coordinates.txt"

# Identifichiamo anche il prefisso comune (es. run=20250826-035846)
# Prendiamo il primo file disponibile per estrarre la parte prima del primo "_"
first_file=$(ls "$ORIGINALS_DIR"/*_data=diff.tif | head -n 1)
common_prefix=$(basename "$first_file" | cut -d'_' -f1)

for f in "$ORIGINALS_DIR"/*_data=diff.tif; do
    filename=$(basename "$f")
    val_T=$(echo "$filename" | sed -n 's/.*T=\([0-9.]*\).*/\1/p')
    if (( $(echo "$val_T > $max_T" | bc -l) )); then
        max_T=$val_T
        reference_file=$f
    fi
done

echo "Max T: $max_T | Prefix: $common_prefix"

# --- FASE 2: GENERAZIONE COORDINATE DI RIFERIMENTO ---
ref_basename=$(basename "$reference_file" .tif)
ref_processed="$PROCESSED_DIR/${ref_basename}_processed.tif"

if [ ! -f "$ref_processed" ]; then
    cd "$BASE_DIR/manipulate_images/"
    python process_image.py --input "$reference_file" --process remove_column_bias remove_hot_pixels remove_cold_pixels --output "$ref_processed"
fi

cd "$BASE_DIR/find_centers/"
python find_centers.py --input "$ref_processed" --convolution --coordinates_root "$reference_coords"

# --- FASE 3: LOOP DI PROCESSAMENTO SU TUTTI I FILE ---
cd "$BASE_DIR"
for input_path in "$ORIGINALS_DIR"/*_data=diff.tif; do
    filename=$(basename "$input_path")
    basename="${filename%.tif}"
    error_file="${input_path%diff.tif}diffe.tif"
    
    if [ ! -f "$error_file" ]; then continue; fi

    val_T=$(echo "$filename" | sed -n 's/.*T=\([0-9.]*\).*/\1/p')
    val_v=$(echo "$filename" | sed -n 's/.*v=\([0-9.]*\).*/\1/p')

    # Processamento (saltiamo se già fatto per il ref)
    if [ ! -f "$PROCESSED_DIR/${basename}_processed.tif" ]; then
        cd "$BASE_DIR/manipulate_images/"
        python process_image.py --input "$input_path" --process remove_column_bias remove_hot_pixels remove_cold_pixels --output "$PROCESSED_DIR/${basename}_processed.tif"
    fi
    
    cd "$BASE_DIR/manipulate_images/"
    python tif2th2.py --input "$PROCESSED_DIR/${basename}_processed.tif" --error "$error_file" --output "$ROOT_DIR/${basename}_processed_th2f.root"

    cd "$BASE_DIR/spot_luminosity/"
    root -l -q "spot_luminosity_sum.C(\"$ROOT_DIR/${basename}_processed_th2f.root\",\"$reference_coords\")"

    if [ -f "luminosity_results.csv" ]; then
        [[ -z "$val_T" || -z "$val_v" ]] && { val_T="0"; val_v="0"; }
        awk -v t="$val_T" -v v="$val_v" 'BEGIN { FS=","; OFS="," } { gsub(/[[:space:]]+$/, "", $0); if ($0 != "") { if (NR == 1) print $0, "T", "v"; else print $0, t, v; } }' luminosity_results.csv > "$LUM_DIR/luminosity_${basename}_processed.csv"
        rm -f luminosity_results.csv
    fi
done

# --- FASE 4: UNIONE, ASSEGNAZIONE SPOT E ORDINAMENTO ---
echo "Merging and sorting results..."
cd "$LUM_DIR"

python3 <<EOF
import pandas as pd
import glob
import os

# 1. Unisce tutti i file csv creati
files = glob.glob("luminosity_*_processed.csv")
df_list = []

for f in files:
    # Leggiamo il file (saltiamo l'intestazione se il dataframe è già iniziato)
    temp_df = pd.read_csv(f)
    df_list.append(temp_df)

full_df = pd.concat(df_list, ignore_index=True)

# 2. Assegna l'ID 'spot' basato sulle coordinate X e Y
# Assumiamo che le colonne si chiamino 'x' e 'y'
# Creiamo una colonna 'coord_group' per raggruppare coppie x,y identiche
full_df['spot'] = full_df.groupby(['x', 'y']).ngroup()

# 3. Ordina per spot e poi per temperatura (T)
# Assicuriamoci che T sia numerico per l'ordinamento
full_df['T'] = pd.to_numeric(full_df['T'], errors='coerce')
full_df = full_df.sort_values(by=['spot', 'T'])

# 4. Salva il file finale
output_name = "${common_prefix}_total.csv"
full_df.to_csv(output_name, index=False)
print(f"Final file created: {output_name}")
EOF

echo "Done."


# root -l 'lum_vs_T.C("run=20250826-035846_v=7.csv")'