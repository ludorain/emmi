#!/bin/bash

set -euo pipefail

# ============================================================
# INPUT
# ============================================================

# Primo input:
# output del codice precedente.
# Deve già contenere global ID nella colonna "spot".
REFERENCE_CSV="DATA_irradiated/merged_files_A1_T=20/A1_T=20_merged_max_vfin_global_ID.csv"

# Secondo input:
# nuovo CSV da aggiungere.
# Qui la colonna "spot" contiene ancora gli ID locali.
NEW_CSV="DATA_irradiated/annealing_T=100_h=25/A1_T=20_20260625-095817/luminosity/A1_T=20_run=20260625-095817_total.csv"


# ============================================================
# PARAMETRI
# ============================================================

MATCH_RADIUS=40.0

# Se false, il codice si ferma se trova nel NEW_CSV uno spot
# che non corrisponde a nessuno spot globale del REFERENCE_CSV.
# Se true, assegna nuovi global ID agli spot non trovati.
ADD_NEW_SPOTS=true

OUTDIR="DATA_irradiated/merged_files_A1_T=20"

MERGED_OUTPUT="$OUTDIR/merged_max_vfin_global_ID_add_100_25.csv"

mkdir -p "$OUTDIR"


# ============================================================
# CONTROLLO FILE
# ============================================================

if [[ ! -f "$REFERENCE_CSV" ]]; then
    echo "Errore: REFERENCE_CSV non trovato:"
    echo "$REFERENCE_CSV"
    exit 1
fi

if [[ ! -f "$NEW_CSV" ]]; then
    echo "Errore: NEW_CSV non trovato:"
    echo "$NEW_CSV"
    exit 1
fi


# ============================================================
# PYTHON
# ============================================================

python3 - "$MATCH_RADIUS" "$OUTDIR" "$REFERENCE_CSV" "$NEW_CSV" "$MERGED_OUTPUT" "$ADD_NEW_SPOTS" <<'PYCODE'

import sys
from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================
# INPUT DA BASH
# ============================================================

match_radius = float(sys.argv[1])
outdir = Path(sys.argv[2])
reference_path = Path(sys.argv[3])
new_path = Path(sys.argv[4])
merged_output = Path(sys.argv[5])
add_new_spots = sys.argv[6].lower() == "true"

outdir.mkdir(parents=True, exist_ok=True)


# ============================================================
# FUNZIONI UTILI
# ============================================================

def extract_phase(path: Path) -> str:
    """
    Estrae il nome della phase dal path.

    Esempio:
    DATA_irradiated/annealing_T=75_h=5/A1_T=20_run=.../luminosity/file.csv

    restituisce:
    annealing_T=75_h=5
    """
    parts = path.parts

    if "DATA_irradiated" in parts:
        idx = parts.index("DATA_irradiated")
        if idx + 1 < len(parts):
            return parts[idx + 1]

    return path.parent.parent.name


def check_columns(df: pd.DataFrame, path: Path):
    required = ["spot", "x", "y", "luminosity", "error", "T", "v", "v_fin"]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            f"\nIl file {path} non contiene le colonne richieste: {missing}\n"
            f"Colonne trovate: {list(df.columns)}"
        )


def detected_to_bool(series: pd.Series) -> pd.Series:
    """
    Converte la colonna detected in booleani,
    anche se viene letta come stringa.
    """
    if series.dtype == bool:
        return series

    return (
        series.astype(str)
              .str.strip()
              .str.lower()
              .isin(["true", "1", "yes", "y"])
    )


def mean_coordinates_per_spot(df: pd.DataFrame) -> pd.DataFrame:
    """
    Riduce le righe dello stesso spot a una sola coordinata media.
    Questa media serve solo per il matching spaziale.
    """
    tmp = df[["spot", "x", "y"]].copy()

    tmp["spot"] = pd.to_numeric(tmp["spot"]).astype(int)
    tmp["x"] = pd.to_numeric(tmp["x"])
    tmp["y"] = pd.to_numeric(tmp["y"])

    return (
        tmp.groupby("spot", as_index=False)
           .agg(
               x_mean=("x", "mean"),
               y_mean=("y", "mean")
           )
    )


def make_new_global_output_name(path: Path) -> Path:
    """
    Crea il nome del CSV completo del NEW_CSV con global ID.

    Esempi:
    file_total.csv -> file_total_global_ID.csv
    file.csv       -> file_total_global_ID.csv
    """
    stem = path.stem

    if stem.endswith("_total"):
        new_name = f"{stem}_global_ID{path.suffix}"
    else:
        new_name = f"{stem}_total_global_ID{path.suffix}"

    return outdir / new_name


def build_catalog_from_reference(reference_df: pd.DataFrame) -> pd.DataFrame:
    """
    Costruisce il catalogo globale a partire dal primo CSV.

    IMPORTANTE:
    La colonna spot del primo CSV è già il global ID.
    """
    df = reference_df.copy()

    # Se c'è detected, uso preferibilmente solo righe realmente viste.
    if "detected" in df.columns:
        detected_mask = detected_to_bool(df["detected"])

        if detected_mask.any():
            df = df[detected_mask].copy()

    ref_spots = mean_coordinates_per_spot(df)

    catalog = ref_spots.rename(columns={
        "spot": "global_spot_id",
        "x_mean": "x_ref",
        "y_mean": "y_ref"
    })

    catalog["global_spot_id"] = catalog["global_spot_id"].astype(int)

    return catalog[["global_spot_id", "x_ref", "y_ref"]]


def greedy_match_to_catalog(local_spots: pd.DataFrame,
                            catalog: pd.DataFrame,
                            radius: float):
    """
    Associa gli spot locali del NEW_CSV agli spot globali del REFERENCE_CSV.

    Ritorna:
    - mapping: local_spot -> global_spot_id
    - distances: local_spot -> distanza dal global spot assegnato
    """

    pairs = []

    for _, loc in local_spots.iterrows():
        local_spot = int(loc["spot"])
        x = loc["x_mean"]
        y = loc["y_mean"]

        for _, ref in catalog.iterrows():
            gid = int(ref["global_spot_id"])
            xr = ref["x_ref"]
            yr = ref["y_ref"]

            d = np.sqrt((x - xr)**2 + (y - yr)**2)

            if d <= radius:
                pairs.append((d, local_spot, gid))

    pairs.sort(key=lambda t: t[0])

    assigned_local = set()
    assigned_global = set()

    mapping = {}
    distances = {}

    for d, local_spot, gid in pairs:
        if local_spot in assigned_local:
            continue

        if gid in assigned_global:
            continue

        mapping[local_spot] = gid
        distances[local_spot] = d

        assigned_local.add(local_spot)
        assigned_global.add(gid)

    return mapping, distances


def print_unmatched_spots(unmatched, local_spots, catalog):
    """
    Stampa per ogni spot non associato il global spot più vicino.
    """
    print("")
    print("Spot del NEW_CSV non associati al riferimento:")
    print("")

    for local_spot in unmatched:
        row = local_spots[local_spots["spot"] == local_spot].iloc[0]
        x = row["x_mean"]
        y = row["y_mean"]

        best_gid = None
        best_dist = None

        for _, ref in catalog.iterrows():
            gid = int(ref["global_spot_id"])
            d = np.sqrt((x - ref["x_ref"])**2 + (y - ref["y_ref"])**2)

            if best_dist is None or d < best_dist:
                best_dist = d
                best_gid = gid

        print(
            f"local_spot={local_spot}, "
            f"x_mean={x:.2f}, y_mean={y:.2f}, "
            f"global più vicino={best_gid}, "
            f"distanza={best_dist:.2f} pixel"
        )

    print("")


def make_reference_templates(reference_df: pd.DataFrame) -> pd.DataFrame:
    """
    Crea una riga template per ogni global spot.
    Serve per aggiungere nel MERGED_OUTPUT righe con luminosity=0,
    error=0 e detected=False quando uno spot non è presente nel NEW_CSV.
    """
    df = reference_df.copy()

    if "detected" in df.columns:
        detected_mask = detected_to_bool(df["detected"])

        if detected_mask.any():
            df = df[detected_mask].copy()

    df["spot"] = pd.to_numeric(df["spot"]).astype(int)
    df["v_fin"] = pd.to_numeric(df["v_fin"])

    idx = df.groupby("spot")["v_fin"].idxmax()

    templates = df.loc[idx].copy()
    templates = templates.set_index("spot")

    return templates


# ============================================================
# LETTURA CSV
# ============================================================

reference_df_original = pd.read_csv(reference_path)
new_df_original = pd.read_csv(new_path)

check_columns(reference_df_original, reference_path)
check_columns(new_df_original, new_path)

new_phase = extract_phase(new_path)

print("CSV di riferimento:")
print(reference_path)
print("")
print("Nuovo CSV da unire:")
print(new_path)
print("")
print(f"Phase del nuovo CSV: {new_phase}")
print(f"Raggio di matching: {match_radius} pixel")


# ============================================================
# COSTRUZIONE CATALOGO DAL PRIMO CSV
# ============================================================

catalog = build_catalog_from_reference(reference_df_original)

print("")
print("Catalogo globale costruito dal primo CSV:")
print(catalog.to_string(index=False))


# ============================================================
# MATCHING DEL NEW_CSV SUL CATALOGO GLOBALE
# ============================================================

local_spots = mean_coordinates_per_spot(new_df_original)

mapping, distances = greedy_match_to_catalog(
    local_spots=local_spots,
    catalog=catalog,
    radius=match_radius
)

all_local_spots = sorted(local_spots["spot"].astype(int).unique())
unmatched = [s for s in all_local_spots if s not in mapping]

if unmatched:
    if not add_new_spots:
        print_unmatched_spots(unmatched, local_spots, catalog)

        raise RuntimeError(
            "Alcuni spot del NEW_CSV non sono stati trovati nel REFERENCE_CSV. "
            "Aumenta MATCH_RADIUS oppure imposta ADD_NEW_SPOTS=true."
        )

    else:
        next_global_id = int(catalog["global_spot_id"].max()) + 1

        for local_spot in unmatched:
            row = local_spots[local_spots["spot"] == local_spot].iloc[0]

            gid = next_global_id
            next_global_id += 1

            mapping[local_spot] = gid
            distances[local_spot] = np.nan

            new_catalog_row = pd.DataFrame([{
                "global_spot_id": gid,
                "x_ref": row["x_mean"],
                "y_ref": row["y_mean"]
            }])

            catalog = pd.concat([catalog, new_catalog_row], ignore_index=True)

            print(
                f"ATTENZIONE: local_spot={local_spot} non era presente "
                f"nel riferimento. Assegnato nuovo global ID = {gid}"
            )


# ============================================================
# OUTPUT 1:
# CSV COMPLETO DEL NEW_CSV CON GLOBAL ID
# ============================================================

new_df_global_full = new_df_original.copy()

# Sostituisce solo la colonna spot.
# Tutte le altre colonne del NEW_CSV restano invariate.
new_df_global_full["spot"] = (
    new_df_global_full["spot"]
    .apply(lambda s: mapping[int(s)])
    .astype(int)
)

# Ordina il CSV completo secondo global ID
sort_cols_full = ["spot"]

if "v_fin" in new_df_global_full.columns:
    sort_cols_full.append("v_fin")
if "T" in new_df_global_full.columns:
    sort_cols_full.append("T")
if "v" in new_df_global_full.columns:
    sort_cols_full.append("v")

new_df_global_full = (
    new_df_global_full
    .sort_values(sort_cols_full)
    .reset_index(drop=True)
)

new_global_output = make_new_global_output_name(new_path)
new_df_global_full.to_csv(new_global_output, index=False)

print("")
print(f"Creato CSV completo con global ID: {new_global_output}")


# ============================================================
# OUTPUT 2:
# MERGED_OUTPUT CON SOLO v_fin MASSIMO DEL NEW_CSV
# ============================================================

new_df_for_summary = new_df_global_full.copy()

new_df_for_summary["spot"] = pd.to_numeric(new_df_for_summary["spot"]).astype(int)
new_df_for_summary["v_fin"] = pd.to_numeric(new_df_for_summary["v_fin"])

# Una sola riga per ogni spot del NEW_CSV: quella a v_fin massimo
idx_new_max = new_df_for_summary.groupby("spot")["v_fin"].idxmax()
new_max_rows = new_df_for_summary.loc[idx_new_max].copy()
new_max_rows = new_max_rows.set_index("spot")

present_spots_new = set(new_max_rows.index.astype(int))

# Tutti gli spot globali presenti nel catalogo di riferimento
all_global_spots = sorted(catalog["global_spot_id"].astype(int).unique())

# Template dal primo CSV, per costruire righe detected=False
reference_templates = make_reference_templates(reference_df_original)

summary_rows_new = []

for gid in all_global_spots:

    if gid in present_spots_new:
        row = new_max_rows.loc[gid].copy()
        row["spot"] = gid

        # Aggiungo/aggiorno phase e detected per il file merged
        row["phase"] = new_phase
        row["detected"] = True

        summary_rows_new.append(row)

    else:
        # Spot presente nel riferimento ma assente nel NEW_CSV:
        # lo aggiungo come non rivelato, come nel codice precedente.
        if gid not in reference_templates.index:
            continue

        row = reference_templates.loc[gid].copy()

        row["spot"] = gid
        row["luminosity"] = 0
        row["error"] = 0
        row["phase"] = new_phase
        row["detected"] = False

        summary_rows_new.append(row)


new_summary = pd.DataFrame(summary_rows_new)


# ============================================================
# MERGE TRA REFERENCE_CSV E NUOVE RIGHE A v_fin MASSIMO
# ============================================================

reference_for_merge = reference_df_original.copy()

# Le colonne finali sono quelle del reference,
# più eventuali colonne extra presenti nel nuovo summary.
final_columns = list(reference_for_merge.columns)

for col in new_summary.columns:
    if col not in final_columns:
        final_columns.append(col)

# Se per qualche motivo phase/detected non fossero nel reference,
# le aggiungo comunque al file merged.
for col in ["phase", "detected"]:
    if col not in final_columns:
        final_columns.append(col)

for col in final_columns:
    if col not in reference_for_merge.columns:
        reference_for_merge[col] = pd.NA

    if col not in new_summary.columns:
        new_summary[col] = pd.NA

reference_for_merge = reference_for_merge[final_columns]
new_summary = new_summary[final_columns]

merged = pd.concat(
    [reference_for_merge, new_summary],
    ignore_index=True
)

# Ordina il file merged per indice globale dello spot,
# come nel codice originale.
sort_cols_merged = ["spot"]

if "phase" in merged.columns:
    sort_cols_merged.append("phase")

merged["spot"] = pd.to_numeric(merged["spot"]).astype(int)

merged = (
    merged
    .sort_values(sort_cols_merged)
    .reset_index(drop=True)
)

merged.to_csv(merged_output, index=False)

print(f"Creato MERGED_OUTPUT: {merged_output}")


# ============================================================
# REPORT FINALE
# ============================================================

print("")
print("============================================")
print("Merge completato.")
print("")
print("Output creati:")
print(f"1) CSV completo NEW_CSV con global ID:")
print(f"   {new_global_output}")
print("")
print(f"2) CSV merged con solo v_fin massimo del NEW_CSV:")
print(f"   {merged_output}")
print("")
print("Mapping applicato al NEW_CSV:")

for local_spot in sorted(mapping.keys()):
    d = distances[local_spot]

    if np.isnan(d):
        d_str = "nuovo spot"
    else:
        d_str = f"{d:.2f} pixel"

    print(
        f"local_spot {local_spot}  ->  "
        f"global_spot_id {mapping[local_spot]}  "
        f"(d = {d_str})"
    )

print("============================================")

PYCODE