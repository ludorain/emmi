#!/bin/bash
#situazione attuale: stiamo lavorando su uno stesso silicio, a cui variamo la temperatura e l'overvoltage
#l'idea è che il numero di difetti resti lo stesso
#l'approccio è quindi di trovare le coordinate dei centri dalla condizione di maggiore visibilità (T maggiore e V maggiore)
#il file con le coordinate dei centri sarà poi usato per tutte le condizioni di temperatura e overvoltage, in modo da confrontare la luminosità degli stessi spot in condizioni diverse

echo "Starting"

# --- PHASE 1: Image cleanup
#parte di pulizia delle immagini
#al momento non implementato perché ce li abbiamo già
#da fare copia e incolla dal file loop_lum_temp.sh

# --- PHASE 2: Creation of th2f files
#parte di creazione dei th2f
#al momento non implementato perché ce li abbiamo già
#da fare copia e incolla dal file loop_lum_temp.sh


#--- PHASE 3: Find defects coordinates
#trovare le coordinate dei centri SOLO sul file con il maggior numero di spot
#quindi quello a T maggiore e V maggiore
cd find_centers/
python find_centers.py \
--input ../files_v=9/2processed/run=20250825-093810_x=95600_y=50600_z=92600_T=22.0_v=9.0_data=diff_processed.tif \
--convolution \
--coordinates_root ../files_v=9/3coordinates/reference_coordinates_type=S13360-3050_run=20250825_26.txt
cd ..

# --- PHASE 4: LUMINOSITY CALCULATION
BASE_DIR=$(pwd)
COORD_DIR="$BASE_DIR/files_v=9/3coordinates"
reference_coords="$COORD_DIR/reference_coordinates_type=S13360-3050_run=20250825_26.txt"

# --- PHASE 4.1: Luminosity calculation for constant v (=9V) and variable T ---
ORIGINALS_DIR1="$BASE_DIR/files_v=9/1originals"
PROCESSED_DIR1="$BASE_DIR/files_v=9/2processed"
ROOT_DIR1="$BASE_DIR/files_v=9/4th2f"
LUM_DIR1="$BASE_DIR/files_v=9/luminosity"
mkdir -p "$PROCESSED_DIR1" "$COORD_DIR" "$ROOT_DIR1" "$LUM_DIR1"

cd "$BASE_DIR/spot_luminosity/"
for input_path in "$ORIGINALS_DIR1"/*_data=diff.tif; do
    filename=$(basename "$input_path")
    basename="${filename%.tif}"
    error_file="${input_path%diff.tif}diffe.tif"
    
    if [ ! -f "$error_file" ]; then continue; fi

    val_T=$(echo "$filename" | sed -n 's/.*T=\([0-9.]*\).*/\1/p')
    val_v=9.0

    root -l -q "spot_luminosity_sum.C(\"$ROOT_DIR1/${basename}_processed_th2f.root\",\"$reference_coords\")"

    if [ -f "luminosity_results.csv" ]; then
        [[ -z "$val_T" || -z "$val_v" ]] && { val_T="0"; val_v="0"; }
        awk -v t="$val_T" -v v="$val_v" 'BEGIN { FS=","; OFS="," } { gsub(/[[:space:]]+$/, "", $0); if ($0 != "") { if (NR == 1) print $0, "T", "v"; else print $0, t, v; } }' luminosity_results.csv > "$LUM_DIR1/luminosity_${basename}_processed.csv"
        rm -f luminosity_results.csv
    fi
done
cd "$BASE_DIR"

# --- PHASE 4.2: Luminosity calculation for constant v (=7V) and variable T ---
ORIGINALS_DIR2="$BASE_DIR/files_v=7/1originals"
PROCESSED_DIR2="$BASE_DIR/files_v=7/2processed"
ROOT_DIR2="$BASE_DIR/files_v=7/4th2f"
LUM_DIR2="$BASE_DIR/files_v=7/luminosity"
mkdir -p "$PROCESSED_DIR2" "$COORD_DIR" "$ROOT_DIR2" "$LUM_DIR2"

cd "$BASE_DIR/spot_luminosity/"
for input_path in "$ORIGINALS_DIR2"/*_data=diff.tif; do
    filename=$(basename "$input_path")
    basename="${filename%.tif}"
    error_file="${input_path%diff.tif}diffe.tif"
    
    if [ ! -f "$error_file" ]; then continue; fi

    val_T=$(echo "$filename" | sed -n 's/.*T=\([0-9.]*\).*/\1/p')
    val_v=7.0

    root -l -q "spot_luminosity_sum.C(\"$ROOT_DIR2/${basename}_processed_th2f.root\",\"$reference_coords\")"

    if [ -f "luminosity_results.csv" ]; then
        [[ -z "$val_T" || -z "$val_v" ]] && { val_T="0"; val_v="0"; }
        awk -v t="$val_T" -v v="$val_v" 'BEGIN { FS=","; OFS="," } { gsub(/[[:space:]]+$/, "", $0); if ($0 != "") { if (NR == 1) print $0, "T", "v"; else print $0, t, v; } }' luminosity_results.csv > "$LUM_DIR2/luminosity_${basename}_processed.csv"
        rm -f luminosity_results.csv
    fi
done
cd "$BASE_DIR"

# --- PHASE 4.3: Luminosity calculation for constant T (=20.0°C) and variable overvoltage ---
ORIGINALS_DIR3="$BASE_DIR/files_T=20/1originals"
PROCESSED_DIR3="$BASE_DIR/files_T=20/2processed"
ROOT_DIR3="$BASE_DIR/files_T=20/4th2f"
LUM_DIR3="$BASE_DIR/files_T=20/luminosity"
mkdir -p "$PROCESSED_DIR3" "$COORD_DIR" "$ROOT_DIR3" "$LUM_DIR3"

cd "$BASE_DIR/spot_luminosity/"
for input_path in "$ORIGINALS_DIR3"/*_data=diff.tif; do
    filename=$(basename "$input_path")
    basename="${filename%.tif}"
    error_file="${input_path%diff.tif}diffe.tif"
    
    if [ ! -f "$error_file" ]; then continue; fi

    val_T=20.0
    # Estraiamo il valore e lasciamo che sia awk a gestire la sottrazione decimale
    v_raw=$(echo "$filename" | sed -n 's/.*v=\([0-9.]*\).*/\1/p')
    val_v=$(awk -v v="$v_raw" 'BEGIN { print v - 51.3 }')

    root -l -q "spot_luminosity_sum.C(\"$ROOT_DIR3/${basename}_processed_th2f.root\",\"$reference_coords\")"

    if [ -f "luminosity_results.csv" ]; then
        [[ -z "$val_T" || -z "$val_v" ]] && { val_T="0"; val_v="0"; }
        
        awk -v t="$val_T" -v v="$val_v" 'BEGIN { FS=","; OFS="," } 
        { 
            sub(/\r$/, "", $0); # Rimuove il carriage return (causa della 7° colonna)
            gsub(/[[:space:]]+$/, "", $0); 
            if ($0 != "") { 
                if (NR == 1) 
                    print $0, "T", "v"; 
                else 
                    printf "%s,%g,%g\n", $0, t, v; 
            } 
        }' luminosity_results.csv > "$LUM_DIR3/luminosity_${basename}_processed.csv"
        
        rm -f luminosity_results.csv
    fi
done


# A questo punto abbiamo i file csv con le luminosità degli spot
# per tutte le condizioni di lavoro, indicizzati solo alle coordinate degli spot
#

# --- PHASE 5: Merge csv in folders
echo "Merging CSV files..."

# Definiamo le cartelle di output della fase precedente
LUM_DIRS=("$LUM_DIR1" "$LUM_DIR2" "$LUM_DIR3")

for dir in "${LUM_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        output_file="$dir/merged_luminosity.csv"
        echo "Processing directory: $dir"
        
        # Rimuove il file di merge se esiste già per evitare append infiniti
        rm -f "$output_file"

        # Prende tutti i file csv della cartella (escludendo il file di merge stesso)
        files=($(ls "$dir"/luminosity_*_processed.csv 2>/dev/null))
        
        if [ ${#files[@]} -eq 0 ]; then
            echo "No CSV files found in $dir"
            continue
        fi

        # Per il primo file, prendiamo tutto (header incluso)
        cat "${files[0]}" > "$output_file"

        # Per i file successivi (dal secondo in poi), saltiamo l'header (NR > 1)
        for (( i=1; i<${#files[@]}; i++ )); do
            awk 'NR > 1' "${files[$i]}" >> "$output_file"
        done
        
        echo "Created: $output_file"
    fi
done

echo "All merges completed."

# --- PHASE 6: Global merge and spot indexing (v=7 and v=9)
echo "Starting global merge and spot indexing with Python..."

python3 <<EOF
import pandas as pd
import os

# Definiamo i percorsi di input basati sulle variabili Bash
dir_v7 = "$LUM_DIR2/merged_luminosity.csv"
dir_v9 = "$LUM_DIR1/merged_luminosity.csv"

# Percorso di output richiesto
output_dir = "$BASE_DIR/spot_luminosity"
output_filename = "run=202508_2overvoltages_scanTemp.csv"
output_path = os.path.join(output_dir, output_filename)

# Assicuriamoci che la cartella di destinazione esista
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Caricamento e processamento
if os.path.exists(dir_v7) and os.path.exists(dir_v9):
    df7 = pd.read_csv(dir_v7)
    df9 = pd.read_csv(dir_v9)
    
    # Unione dei dataframe
    df_global = pd.concat([df7, df9], ignore_index=True)
    
    # Pulizia nomi colonne
    df_global.columns = df_global.columns.str.strip()
    
    # Creazione variabile 'spot' basata su coordinate univoche
    # ngroup() assegna lo stesso ID agli stessi (x, y)
    df_global['spot'] = df_global.groupby(['x', 'y']).ngroup()
    
    # Riordino colonne: spot come prima colonna
    cols = ['spot'] + [c for c in df_global.columns if c != 'spot']
    df_global = df_global[cols]
    
    # Ordinamento richiesto: spot -> T -> v
    df_global = df_global.sort_values(by=['spot', 'T', 'v'])
    
    # Salvataggio in spot_luminosity/run=202508_2overvoltages_scanTemp.csv
    df_global.to_csv(output_path, index=False)
    print(f"Success! Final merge saved in: {output_path}")
else:
    if not os.path.exists(dir_v7): print(f"Missing file: {dir_v7}")
    if not os.path.exists(dir_v9): print(f"Missing file: {dir_v9}")
EOF


# --- PHASE 7: Merge T=20 and apply consistent spot ID
echo "Starting Phase 7: ScanV merge and consistent indexing..."

python3 <<EOF
import pandas as pd
import os
import glob

# Percorsi
dir_t20 = "$LUM_DIR3"
reference_file = "$BASE_DIR/spot_luminosity/run=202508_2overvoltages_scanTemp.csv"
output_dir_spot = "$BASE_DIR/spot_luminosity"

# 1) Merge dei file nella cartella T=20
all_files = glob.glob(os.path.join(dir_t20, "luminosity_*_processed.csv"))
if not all_files:
    print("Error: No CSV files found in $LUM_DIR3")
    exit()

df_list = [pd.read_csv(f) for f in all_files]
df_scanV = pd.concat(df_list, ignore_index=True)
df_scanV.columns = df_scanV.columns.str.strip()

# 2) Salva il merge grezzo nella cartella T=20
scanV_raw_path = os.path.join(dir_t20, "run=202508_scanV.csv")
df_scanV.to_csv(scanV_raw_path, index=False)
print(f"Merged T=20 file saved: {scanV_raw_path}")

# 3) Creazione copia con variabile 'spot' coerente
if os.path.exists(reference_file):
    # Carichiamo il file della Phase 6 per estrarre la mappa coordinate -> spot
    df_ref = pd.read_csv(reference_file)
    # Creiamo un set univoco di x, y e il relativo spot
    mapping = df_ref[['x', 'y', 'spot']].drop_duplicates()
    
    # Uniamo il mapping al nuovo dataframe del T=20
    # Usiamo 'left' merge per mantenere tutte le righe di scanV
    df_scanV_indexed = pd.merge(df_scanV, mapping, on=['x', 'y'], how='left')
    
    # Riordino colonne per avere 'spot' all'inizio
    cols = ['spot'] + [c for c in df_scanV_indexed.columns if c != 'spot']
    df_scanV_indexed = df_scanV_indexed[cols]
    
    # Ordinamento: spot, poi T, poi v
    df_scanV_indexed = df_scanV_indexed.sort_values(by=['spot', 'T', 'v'])
    
    # Salvataggio nella cartella spot_luminosity
    final_output_path = os.path.join(output_dir_spot, "run=202508_scanV.csv")
    df_scanV_indexed.to_csv(final_output_path, index=False)
    print(f"Consistent indexed file saved: {final_output_path}")
else:
    print(f"Error: Reference file {reference_file} not found. Cannot sync spot IDs.")
EOF