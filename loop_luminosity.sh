#!/bin/bash

BASE_DIR=$(pwd)
ORIGINALS_DIR="$BASE_DIR/files/1originals"
PROCESSED_DIR="$BASE_DIR/files/2processed"
COORD_DIR="$BASE_DIR/files/3coordinates"
ROOT_DIR="$BASE_DIR/files/4th2f"
LUM_DIR="$BASE_DIR/files/luminosity"

mkdir -p "$PROCESSED_DIR" "$COORD_DIR" "$ROOT_DIR" "$LUM_DIR"

for input_path in "$ORIGINALS_DIR"/*_data=diff.tif; do
    
    filename=$(basename "$input_path")
    basename="${filename%.tif}"
    error_file="${input_path%diff.tif}diffe.tif"
    
    if [ ! -f "$error_file" ]; then
        echo "Warning: Error file not found for $filename, skipping."
        continue
    fi

    # --- NUOVO METODO ESTRAZIONE (PIÙ COMPATIBILE) ---
    # Estraiamo T: cerchiamo T=, prendiamo tutto ciò che segue fino al prossimo _
    val_T=$(echo "$filename" | sed -n 's/.*T=\([0-9.]*\).*/\1/p')
    # Estraiamo v: cerchiamo v=, prendiamo tutto ciò che segue fino al prossimo _
    val_v=$(echo "$filename" | sed -n 's/.*v=\([0-9.]*\).*/\1/p')

    echo "Processing: $filename (T=$val_T, v=$val_v)"

    # 1, 2, 3, 4: Esecuzione script (omessi per brevità, rimangono uguali a prima)
    cd "$BASE_DIR/manipulate_images/"
    python process_image.py --input "$input_path" --process remove_column_bias remove_hot_pixels remove_cold_pixels --output "$PROCESSED_DIR/${basename}_processed.tif"
    
    python tif2th2.py --input "$PROCESSED_DIR/${basename}_processed.tif" --error "$error_file" --output "$ROOT_DIR/${basename}_processed_th2f.root"

    cd "$BASE_DIR/find_centers/"
    python find_centers.py --input "$PROCESSED_DIR/${basename}_processed.tif" --convolution --coordinates_root "$COORD_DIR/${basename}_processed_coords.txt"

    cd "$BASE_DIR/spot_luminosity/"
    root -l -q "spot_luminosity_radius.C(\"$ROOT_DIR/${basename}_processed_th2f.root\",\"$COORD_DIR/${basename}_processed_coords.txt\")"

    # --- MODIFICA CSV CON INTESTAZIONE ---
    if [ -f "luminosity_results.csv" ]; then
        # Verifichiamo se le variabili non sono vuote
        if [[ -z "$val_T" || -z "$val_v" ]]; then
            val_T="NaN"
            val_v="NaN"
        fi

        # NR == 1 è la prima riga (header)
        # NR > 1 sono le righe con i dati
        awk -v t="$val_T" -v v="$val_v" '
        BEGIN { FS=","; OFS="," }
        { 
            if (NF > 0 && $0 !~ /^[[:space:]]*$/) {
                if (NR == 1) {
                    print $0, "T", "v"
                } else {
                    print $0, t, v
                }
            }
        }' luminosity_results.csv > "temp_results.csv"

        mv "temp_results.csv" "$LUM_DIR/luminosity_${basename}_processed.csv"
        rm -f luminosity_results.csv
        echo "Successfully saved to $LUM_DIR with headers."
    else
        echo "Error: luminosity_results.csv not found."
    fi

done