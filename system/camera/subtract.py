#! /usr/bin/env python

import argparse
import tifffile
import numpy as np

def subtract_tifs(input1, input2, output):
    # Load both images
    image1 = tifffile.imread(input1).astype(np.float32)
    image2 = tifffile.imread(input2).astype(np.float32)

    # Check if dimensions match
    if image1.shape != image2.shape:
        raise ValueError("Error: Input images must have the same dimensions.")

    # Subtract images
    result = image1 - image2

    # Save the result
    tifffile.imwrite(output, result.astype(np.float32))  # Convert to int16 to allow negative values
    print(f"Subtracted image saved to {output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subtract two TIF images.")
    parser.add_argument("--input1", required=True, help="Path to the first input TIF file.")
    parser.add_argument("--input2", required=True, help="Path to the second input TIF file.")
    parser.add_argument("--output", required=True, help="Path to save the output TIF file.")

    args = parser.parse_args()
    subtract_tifs(args.input1, args.input2, args.output)
