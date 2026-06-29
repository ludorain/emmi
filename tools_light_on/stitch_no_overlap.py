#! /usr/bin/env python

import argparse
import os
import re
import numpy as np
import tifffile
import matplotlib.pyplot as plt


def extract_coords(filename):
    match = re.search(r"x=(-?\d+(?:\.\d+)?)\.y=(-?\d+(?:\.\d+)?)", filename)
    if not match:
        raise ValueError(f"Coordinate non trovate in {filename}")
    
    x = float(match.group(1))
    y = float(match.group(2))
    return x, y


def main():
    parser = argparse.ArgumentParser(description="Stitch immagini con orientazione fisica corretta")
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output")
    parser.add_argument("--display", action="store_true")
    args = parser.parse_args()

    files = [f for f in os.listdir(args.input_dir) if f.lower().endswith(".tif")]

    if len(files) == 0:
        raise ValueError("Nessun file .tif trovato")

    images = {}
    xs, ys = [], []

    # Lettura immagini
    for f in files:
        x, y = extract_coords(f)
        img = tifffile.imread(os.path.join(args.input_dir, f))

        if img.ndim == 3:
            img = img.mean(axis=2)

        images[(x, y)] = img
        xs.append(x)
        ys.append(y)

    #  ORDINAMENTO CORRETTO (nuova convenzione)
    xs_sorted = sorted(set(xs), reverse=True)  # x cresce → sinistra
    ys_sorted = sorted(set(ys), reverse=True)  # y cresce → alto

    nx = len(xs_sorted)
    ny = len(ys_sorted)

    print(f"Grid: {nx} x {ny}")

    #  Dimensione immagini
    sample = next(iter(images.values()))
    h, w = sample.shape

    stitched = np.zeros((ny * h, nx * w), dtype=sample.dtype)

    #  Riempimento
    for (x, y), img in images.items():

        col = xs_sorted.index(x)  # x → colonna
        row = ys_sorted.index(y)  # y → riga

        y0 = row * h
        y1 = y0 + h

        x0 = col * w
        x1 = x0 + w

        stitched[y0:y1, x0:x1] = img

     # Laboratory view (180° rotation)
    visuale_lab = np.rot90(stitched, 2)


    # Salvataggio
    if args.output:
        tifffile.imwrite(args.output, stitched)
        

    if args.display:
        # First canvas
        plt.figure()
        plt.imshow(stitched)
        plt.title("Camera frame of reference")
        plt.axis("off")

        # Second canvas
        plt.figure()
        plt.imshow(visuale_lab)
        plt.title("Laboratory frame of reference")
        plt.axis("off")

        plt.show()

if __name__ == "__main__":
    main()