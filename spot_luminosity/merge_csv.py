#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np


REQUIRED_COLUMNS = {"x", "y", "luminosity", "error", "T", "v", "spot"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge di due CSV assegnando nuovi spot in base alle coordinate x-y."
    )
    parser.add_argument("csv1", help="Primo file CSV")
    parser.add_argument("csv2", help="Secondo file CSV")
    parser.add_argument("-o", "--output", default="merged_output.csv", help="File CSV di output")
    parser.add_argument(
        "--tol",
        type=float,
        default=4.0,
        help="Tolleranza massima sulle coordinate x e y per considerare due spot uguali (default: 4.0)"
    )
    return parser.parse_args()


def check_columns(df, filename):
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"Nel file '{filename}' mancano le colonne richieste: {sorted(missing)}"
        )


def detect_file_voltage(df, filename):
    """
    Assume che in ciascun file il valore di v sia costante (come da tua descrizione).
    """
    unique_v = sorted(df["v"].dropna().unique())
    if len(unique_v) != 1:
        raise ValueError(
            f"Il file '{filename}' dovrebbe avere un solo valore di v, ma ne ha {unique_v}"
        )
    return float(unique_v[0])


def coords_match(x1, y1, x2, y2, tol=4.0):
    """
    Due coordinate sono considerate uguali se |x1-x2| <= tol e |y1-y2| <= tol.
    """
    return abs(x1 - x2) <= tol and abs(y1 - y2) <= tol


def assign_spot_ids(df_all, tol=4.0):
    """
    Assegna nuovi spot ignorando i valori originali.
    Strategia:
    - si scorrono le righe
    - ogni riga viene assegnata al primo gruppo con centro compatibile entro tol
    - altrimenti si crea un nuovo gruppo

    Il centro del gruppo viene aggiornato come media delle coordinate assegnate.
    """
    groups = []
    assigned_spots = []

    for _, row in df_all.iterrows():
        x, y = float(row["x"]), float(row["y"])
        assigned = False

        for group in groups:
            if coords_match(x, y, group["x_mean"], group["y_mean"], tol):
                group["count"] += 1
                group["x_mean"] = (group["x_mean"] * (group["count"] - 1) + x) / group["count"]
                group["y_mean"] = (group["y_mean"] * (group["count"] - 1) + y) / group["count"]
                assigned_spots.append(group["spot"])
                assigned = True
                break

        if not assigned:
            new_spot = len(groups)
            groups.append({
                "spot": new_spot,
                "x_mean": x,
                "y_mean": y,
                "count": 1
            })
            assigned_spots.append(new_spot)

    df_all = df_all.copy()
    df_all["spot"] = assigned_spots
    return df_all


def build_spot_reference(df_all):
    """
    Per ogni nuovo spot calcola coordinate rappresentative (medie).
    """
    ref = (
        df_all.groupby("spot", as_index=False)
        .agg({
            "x": "mean",
            "y": "mean"
        })
        .rename(columns={"x": "x_ref", "y": "y_ref"})
    )
    return ref


def complete_missing_rows(df, ref_df, file_v):
    """
    Per ogni spot e per ogni T globale, se nel file manca una riga,
    la crea con luminosity=0 e error=0.
    """
    all_T = sorted(ref_df["T"].unique())
    all_spots = sorted(ref_df["spot"].unique())

    existing = set(zip(df["spot"], df["T"]))
    missing_rows = []

    for spot in all_spots:
        spot_ref = ref_df[ref_df["spot"] == spot].iloc[0]
        x_ref = spot_ref["x_ref"]
        y_ref = spot_ref["y_ref"]

        for T in all_T:
            if (spot, T) not in existing:
                missing_rows.append({
                    "x": x_ref,
                    "y": y_ref,
                    "luminosity": 0.0,
                    "error": 0.0,
                    "T": T,
                    "v": file_v,
                    "spot": spot
                })

    if missing_rows:
        df_missing = pd.DataFrame(missing_rows)
        df = pd.concat([df, df_missing], ignore_index=True)

    return df


def main():
    args = parse_args()

    df1 = pd.read_csv(args.csv1)
    df2 = pd.read_csv(args.csv2)

    check_columns(df1, args.csv1)
    check_columns(df2, args.csv2)

    v1 = detect_file_voltage(df1, args.csv1)
    v2 = detect_file_voltage(df2, args.csv2)

    df1 = df1.copy()
    df2 = df2.copy()

    df1["source_file"] = 1
    df2["source_file"] = 2

    # Ignora spot originali e riassegna dopo il merge
    df_all = pd.concat([df1, df2], ignore_index=True)

    # Riassegna spot globali usando x,y
    df_all = assign_spot_ids(df_all, tol=args.tol)

    # Separa di nuovo i due file
    df1_new = df_all[df_all["source_file"] == 1].copy()
    df2_new = df_all[df_all["source_file"] == 2].copy()

    # Coordinate di riferimento per ogni spot
    spot_ref = build_spot_reference(df_all)

    # Ref con tutti i T globali
    ref_with_T = df_all[["spot", "T"]].drop_duplicates().merge(spot_ref, on="spot", how="left")

    # Completa i buchi in ciascun file
    df1_new = complete_missing_rows(df1_new, ref_with_T, v1)
    df2_new = complete_missing_rows(df2_new, ref_with_T, v2)

    # Merge finale
    merged = pd.concat([df1_new, df2_new], ignore_index=True)

    # Ordine colonne
    merged = merged[["x", "y", "luminosity", "error", "T", "v", "spot"]]

    # Ordinamento richiesto
    merged = merged.sort_values(by=["spot", "T", "v"], ascending=[True, True, True])

    # Reset index
    merged = merged.reset_index(drop=True)

    merged.to_csv(args.output, index=False)
    print(f"File salvato in: {args.output}")


if __name__ == "__main__":
    main()