#! /usr/bin/env python

import argparse
import skimage
import numpy as np

parser = argparse.ArgumentParser(description='Process an EMMI image')
parser.add_argument('--input', type=str, required=True, help='Input TIF filename')
args = parser.parse_args()

img = skimage.io.imread(args.input, as_gray=True)
lap = skimage.filters.laplace(img)

var = lap.var()
N = lap.size
sigma = var * np.sqrt(2/(N-1))

print(var, sigma)
