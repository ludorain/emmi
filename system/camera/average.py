#! /usr/bin/env python

import argparse
import tifffile
import numpy as np
from astropy.stats import sigma_clip

def parse_arguments():
    parser = argparse.ArgumentParser(description="Compute the average image from a stacked TIF file")
    parser.add_argument("--input", required=True, help="Path to the input stacked TIF file")
    parser.add_argument("--output", required=True, help="Path to save the output average TIF file")
    parser.add_argument("--output_error", required=False, help="Path to save the output average error TIF file")
    parser.add_argument("--filter", action="store_true", help="Filter outliers with astropy.stats.sigma_clip")
    return parser.parse_args()

def read_tif(filename):
    images = []
    with tifffile.TiffFile(filename) as tif:
        print(f' --- reading images from {filename}: found {len(tif.pages)} frames')
        for i, page in enumerate(tif.pages):
            image = page.asarray().astype(np.float32)
            images.append(image)
    return np.array(images)
    
def average_tif(images, filter):
    if filter:
        print(f' --- filtering images from outliers')
        images = sigma_clip(images, axis=0)
    average = np.mean(images, axis=0)
    error = np.std(images, axis=0) / np.sqrt(images.shape[0])
    return average.astype(np.float32), error.astype(np.float32)

def write_tif(filename, images):
    tifffile.imwrite(filename, images)
    print(f' --- image saved to {filename}')

if __name__ == "__main__":
    args = parse_arguments()
    images = read_tif(args.input)
    average, error = average_tif(images, args.filter)
    write_tif(args.output, average)
    if args.output_error is not None:
        write_tif(args.output_error, error)
