#!/bin/bash

LUM_DIR="files/luminosity"
OUTPUT_FILE="total_temperature.csv"
TEMP_ALL="all_data_tmp.csv"

echo "Merging and sorting CSV files..."

# 1. Unione dei file
# Prendiamo l'header dal primo file trovato
first_file=$(ls "$LUM_DIR"/*.csv | head -n 1)
head -n 1 "$first_file" > "$TEMP_ALL"

# Accodiamo i dati di tutti i file saltando la prima riga (l'header)
for f in "$LUM_DIR"/*.csv; do
    tail -n +2 "$f" >> "$TEMP_ALL"
done

# 2. Ordinamento
# Supponendo che nel tuo CSV:
# - La colonna "spot" sia la prima ($1) -> opzione -k1,1
# - La colonna "T" sia la penultima (nella versione precedente l'abbiamo messa in coda)
# Se T e v sono le ultime due, e ipotizzando che il CSV originale avesse ad esempio 5 colonne,
# T sarà la colonna 6 e v la colonna 7.
# Nel comando sort, -t',' definisce la virgola come separatore.

# NOTA: Regola i numeri dopo -k in base alla posizione reale delle colonne.
# Se 'spot' è la colonna 1 e 'T' è la colonna 6:
header=$(head -n 1 "$TEMP_ALL")
tail -n +2 "$TEMP_ALL" | sort -t',' -k1,1n -k6,6n > "sorted_body.csv"

# 3. Ricomposizione finale
echo "$header" > "$OUTPUT_FILE"
cat "sorted_body.csv" >> "$OUTPUT_FILE"

# Pulizia file temporanei
rm "$TEMP_ALL" "sorted_body.csv"

echo "Done! File created: $OUTPUT_FILE"