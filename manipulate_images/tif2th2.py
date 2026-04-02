#! /usr/bin/env pyroot

import argparse
import ROOT
import numpy as np
from PIL import Image


def read_tif_as_array(filename):
    img = Image.open(filename)
    arr = np.array(img)

    if arr.ndim != 2:
        raise RuntimeError(f"{filename}: expected a 2D grayscale TIFF, got shape {arr.shape}")

    return arr


def tif_to_th2f(content_file, error_file, histname="h2", histtitle="h2", invert_y=True):
    content = read_tif_as_array(content_file).astype(np.float64)
    if error_file is not None:
        error   = read_tif_as_array(error_file).astype(np.float64)
        if content.shape != error.shape:
            raise RuntimeError(
                f"Shape mismatch: content TIFF has shape {content.shape}, "
                f"error TIFF has shape {error.shape}"
            )
    else:
        error = None

    ny, nx = content.shape

    h2 = ROOT.TH2F(histname, histtitle, nx, 0, nx, ny, 0, ny)
    h2.Sumw2()  # make sure the error storage exists

    for y in range(ny):
        for x in range(nx):
            binx = x + 1
            biny = ny - y if invert_y else y + 1

            val = float(content[y, x])
            err = float(error[y, x]) if error is not None else 0.0
            
            h2.SetBinContent(binx, biny, val)
            h2.SetBinError(binx, biny, err)

    return h2


def main():

    parser = argparse.ArgumentParser(description='TIFF to TH2F converter')
    parser.add_argument('--input', type=str, required=True, help='Input TIF filename of the image')
    parser.add_argument('--error', type=str, required=False, help='Input TIF filename of the errors')
    parser.add_argument('--output', type=str, required=True, help='Output ROOT filename')
    parser.add_argument('--name', type=str, default='h2', required=False, help='TH2F object name')
    parser.add_argument('--title', type=str, default='h2', required=False, help='TH2F object title')
    args = parser.parse_args()

    content_file = args.input
    error_file   = args.error
    output_file  = args.output
    histname     = args.name
    histtitle    = args.title

    h2 = tif_to_th2f(content_file, error_file, histname, histtitle, invert_y=True)

    fout = ROOT.TFile(output_file, "RECREATE")
    h2.Write()
    fout.Close()

    print(f"Wrote histogram '{histname}' to {output_file}")


if __name__ == "__main__":
    main()
