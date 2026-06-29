#Qui prende i file CSV di luminosità dei vari dataset, li confronta e assegna un ID globale agli spot, 
#salvando i file con gli ID globali (che servono per le analisi lum vs T/V)
# e crea un file con UNA RIGA PER OGNI SPOT E PER OGNI PHASE,
# PRENDENDO SOLO v_fin MASSIMO (serve per vedere risultati annealing)


#!/bin/bash

set -euo pipefail

# ============================================================
# INSERISCI QUI I 4 PATH AI CSV
# ============================================================

CSV_1="DATA_irradiated/before_annealing/A1_v=5_run=20260514-062228/luminosity/A1_v=5_run=20260514-062228_total.csv"

CSV_2="DATA_irradiated/annealing_T=75_h=5/A1_v=5_run=20260521-080451/luminosity/A1_v=5_run=20260521-080451_total.csv"

CSV_3="DATA_irradiated/annealing_T=75_h=25/A1_v=5_run=20260531-095134/luminosity/A1_v=5_run=20260531-095134_total.csv"

CSV_4="DATA_irradiated/annealing_T=100_h=5/A1_v=5_run=20260612-001017/luminosity/A1_v=5_run=20260612-001017_total.csv"


# ============================================================
# PARAMETRI
# ============================================================

# Tolleranza spaziale in pixel per dire che due spot sono lo stesso spot.
MATCH_RADIUS=40.0

# Cartella di output
OUTDIR="DATA_irradiated/merged_files_v=const"

mkdir -p "$OUTDIR"


# ============================================================
# CONTROLLO FILE
# ============================================================

CSV_PATHS=("$CSV_1" "$CSV_2" "$CSV_3" "$CSV_4")

for file in "${CSV_PATHS[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "Errore: file non trovato:"
        echo "$file"
        exit 1
    fi
done


# ============================================================
# PYTHON
# ============================================================

python3 - "$MATCH_RADIUS" "$OUTDIR" "${CSV_PATHS[@]}" <<'PYCODE'

import sys
from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================
# INPUT DA BASH
# ============================================================

match_radius = float(sys.argv[1])
outdir = Path(sys.argv[2])
csv_paths = [Path(p) for p in sys.argv[3:]]

outdir.mkdir(parents=True, exist_ok=True)


# ============================================================
# FUNZIONI UTILI
# ============================================================

def extract_phase(path: Path) -> str:
    """
    Estrae il nome della phase dal path.

    Esempio:
    DATA_irradiated/before_annealing/A1_v=5_run=.../luminosity/file.csv

    restituisce:
    before_annealing
    """
    parts = path.parts

    if "DATA_irradiated" in parts:
        idx = parts.index("DATA_irradiated")
        if idx + 1 < len(parts):
            return parts[idx + 1]

    return path.parent.parent.name


def check_columns(df: pd.DataFrame, path: Path):
    """
    Controlla che i CSV abbiano esattamente le colonne necessarie.
    Formato atteso:
    spot,x,y,luminosity,error,T,v
    """
    required = ["spot", "x", "y", "luminosity", "error", "T", "v"]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            f"\nIl file {path} non contiene le colonne richieste: {missing}\n"
            f"Colonne trovate: {list(df.columns)}"
        )


def make_output_name(path: Path) -> Path:
    """
    Crea il nome del file global_ID.
    Esempio:
    file_total.csv -> file_total_global_ID.csv
    """
    return outdir / f"{path.stem}_global_ID{path.suffix}"


def mean_coordinates_per_local_spot(df: pd.DataFrame) -> pd.DataFrame:
    """
    Riduce le righe dello stesso spot locale a una sola coordinata media.
    Questo serve solo per fare il matching.
    I file finali conserveranno comunque tutte le righe originali.
    """
    return (
        df.groupby("spot", as_index=False)
          .agg(
              x_mean=("x", "mean"),
              y_mean=("y", "mean")
          )
    )


def greedy_match_to_catalog(local_spots: pd.DataFrame, catalog: pd.DataFrame, radius: float):
    """
    Associa gli spot locali al catalogo globale usando la distanza spaziale.

    Ritorna:
    - mapping: local_spot -> global_spot_id
    - distances: local_spot -> distanza dal global spot assegnato
    """

    pairs = []

    for _, loc in local_spots.iterrows():
        local_spot = loc["spot"]
        x = loc["x_mean"]
        y = loc["y_mean"]

        for _, ref in catalog.iterrows():
            gid = ref["global_spot_id"]
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


# ============================================================
# LETTURA FILE
# ============================================================

datasets = []

for path in csv_paths:
    df = pd.read_csv(path)
    check_columns(df, path)

    # Forza tutta la colonna v al valore 5
    df["v"] = "5"

    phase = extract_phase(path)

    df["__phase"] = phase
    df["__path"] = str(path)

    datasets.append({
        "path": path,
        "phase": phase,
        "df": df
    })


# ============================================================
# IDENTIFICAZIONE DEL DATASET DI RIFERIMENTO
# ============================================================

reference_dataset = None

for d in datasets:
    if d["phase"] == "before_annealing":
        reference_dataset = d
        break

if reference_dataset is None:
    reference_dataset = datasets[0]
    print(
        "ATTENZIONE: nessun file con phase = before_annealing trovato. "
        "Uso il primo CSV come riferimento."
    )

reference_phase = reference_dataset["phase"]

print(f"Dataset di riferimento: {reference_phase}")
print(f"Raggio di matching: {match_radius} pixel")


# ============================================================
# COSTRUZIONE CATALOGO GLOBALE A PARTIRE DAL BEFORE_ANNEALING
# ============================================================

ref_df = reference_dataset["df"]
ref_local = mean_coordinates_per_local_spot(ref_df)

ref_local = ref_local.sort_values("spot").reset_index(drop=True)
ref_local["global_spot_id"] = np.arange(len(ref_local), dtype=int)

catalog = ref_local.rename(columns={
    "x_mean": "x_ref",
    "y_mean": "y_ref"
})[["global_spot_id", "x_ref", "y_ref"]].copy()

next_global_id = int(catalog["global_spot_id"].max()) + 1 if len(catalog) > 0 else 0


# ============================================================
# ASSEGNAZIONE GLOBAL ID A TUTTI I DATASET
# ============================================================

matched_datasets = []

for d in datasets:
    path = d["path"]
    phase = d["phase"]
    df = d["df"].copy()

    local_spots = mean_coordinates_per_local_spot(df)

    if phase == reference_phase:
        mapping = dict(zip(ref_local["spot"], ref_local["global_spot_id"]))
        distances = {local_spot: 0.0 for local_spot in mapping.keys()}

    else:
        mapping, distances = greedy_match_to_catalog(
            local_spots=local_spots,
            catalog=catalog,
            radius=match_radius
        )

        # Eventuali spot non matched vengono considerati nuovi spot globali.
        # Questo evita di perderli nei file global_ID.
        for _, row in local_spots.iterrows():
            local_spot = row["spot"]

            if local_spot not in mapping:
                gid = next_global_id
                next_global_id += 1

                mapping[local_spot] = gid
                distances[local_spot] = np.nan

                new_row = pd.DataFrame([{
                    "global_spot_id": gid,
                    "x_ref": row["x_mean"],
                    "y_ref": row["y_mean"]
                }])

                catalog = pd.concat([catalog, new_row], ignore_index=True)

                print(
                    f"ATTENZIONE: nuovo spot non presente nel riferimento. "
                    f"phase={phase}, local_spot={local_spot}, "
                    f"assegnato global_spot_id={gid}"
                )

    # Sostituisce la colonna spot con l'ID globale.
    df["spot"] = df["spot"].map(mapping).astype(int)

    # Rimuove le colonne tecniche prima del salvataggio.
    df_to_save = df.drop(columns=["__phase", "__path"], errors="ignore")

    # Forza di nuovo v=5 per sicurezza.
    df_to_save["v"] = "5"

    # Ordina il CSV secondo lo spot globale.
    sort_cols = ["spot"]

    if "T" in df_to_save.columns:
        sort_cols.append("T")
    if "v" in df_to_save.columns:
        sort_cols.append("v")

    df_to_save = df_to_save.sort_values(sort_cols).reset_index(drop=True)

    output_path = make_output_name(path)
    df_to_save.to_csv(output_path, index=False)

    print(f"Creato: {output_path}")

    matched_datasets.append({
        "path": path,
        "phase": phase,
        "df": df_to_save
    })


# ============================================================
# CREAZIONE CSV RIASSUNTIVO
# UNA RIGA PER OGNI SPOT E PER OGNI PHASE,
# PRENDENDO SOLO T MASSIMO
# ============================================================

before_matched = None

for d in matched_datasets:
    if d["phase"] == reference_phase:
        before_matched = d["df"].copy()
        break

if before_matched is None:
    raise RuntimeError("Errore interno: dataset di riferimento non trovato dopo il matching.")


# Per ogni spot globale nel riferimento, prendo la riga con T massimo.
idx_before_max = before_matched.groupby("spot")["T"].idxmax()
before_templates = before_matched.loc[idx_before_max].copy()
before_templates = before_templates.set_index("spot")


# Lista completa degli spot globali.
# Include sia quelli del before_annealing sia eventuali nuovi spot apparsi dopo.
all_global_spots = sorted(catalog["global_spot_id"].astype(int).unique())


summary_rows = []

for d in matched_datasets:
    phase = d["phase"]
    df = d["df"].copy()

    # Forza v=5 anche nel summary.
    df["v"] = "5"

    present_spots = set(df["spot"].astype(int).unique())

    # Riga con T massimo per ogni spot presente in questa phase.
    idx_max = df.groupby("spot")["T"].idxmax()
    max_rows = df.loc[idx_max].copy()
    max_rows = max_rows.set_index("spot")

    for gid in all_global_spots:

        if gid in present_spots:
            row = max_rows.loc[gid].copy()

            # Reinserisce esplicitamente lo spot globale come colonna.
            row["spot"] = gid
            row["v"] = "5"

            row["phase"] = phase
            row["detected"] = True
            summary_rows.append(row)

        else:
            # Caso principale:
            # lo spot era nel before_annealing ma non è presente in questa phase.
            if gid in before_templates.index:
                row = before_templates.loc[gid].copy()

            else:
                # Caso raro:
                # spot nuovo apparso dopo il before_annealing.
                # Non esiste una riga template nel before, quindi uso la prima riga
                # disponibile da un altro dataset.
                template_found = False

                for other in matched_datasets:
                    other_df = other["df"]
                    possible = other_df[other_df["spot"] == gid]

                    if len(possible) > 0:
                        row = possible.loc[possible["T"].idxmax()].copy()
                        template_found = True
                        break

                if not template_found:
                    continue

            row["spot"] = gid
            row["luminosity"] = 0
            row["error"] = 0
            row["v"] = "5"
            row["phase"] = phase
            row["detected"] = False

            summary_rows.append(row)


summary = pd.DataFrame(summary_rows)

# Ordina il file merged per spot globale e phase.
summary = summary.sort_values(["spot", "phase"]).reset_index(drop=True)


# Riordino colonne: originali + phase + detected.
base_columns = [c for c in matched_datasets[0]["df"].columns]
final_columns = base_columns + ["phase", "detected"]

summary = summary[final_columns]

summary_output = outdir / "merged_max_T_global_ID.csv"
summary.to_csv(summary_output, index=False)

print(f"Creato: {summary_output}")


# ============================================================
# REPORT FINALE
# ============================================================

print("")
print("============================================")
print("Output completato.")
print("File creati nella cartella:")
print(outdir)
print("")
print("Catalogo globale:")
print(catalog.to_string(index=False))
print("============================================")

PYCODE