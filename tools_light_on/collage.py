#! /usr/bin/env python

from PIL import Image
import argparse
import math

parser = argparse.ArgumentParser(description='Make a collage of images')
parser.add_argument('--input', type=str, required=True, help='List of image filenames')
parser.add_argument('--output', type=str, required=True, help='Output filename')
parser.add_argument('--columns', type=int, required=True, help='Number of columns')
parser.add_argument('--scale', type=float, default=0.5, required=False, help='Scaling factor')
parser.add_argument('--border', type=int, default=4, required=False, help='Border size in pixels')
args = parser.parse_args()  

cols = args.columns
scale = args.scale
border = args.border

with open(args.input) as f:
    paths = [line.strip() for line in f if line.strip()]

tiles = []

for p in paths:
    im = Image.open(p).convert("RGB")
    w, h = im.size
    new_size = (int(w * scale), int(h * scale))
    im = im.resize(new_size, Image.Resampling.LANCZOS)

    # create tile with white border
    tile = Image.new("RGB", (new_size[0] + 2*border, new_size[1] + 2*border), "white")
    tile.paste(im, (border, border))
    tiles.append(tile)

tw, th = tiles[0].size
rows = math.ceil(len(tiles) / cols)

collage_w = cols * tw
collage_h = rows * th

collage = Image.new("RGB", (collage_w, collage_h), "white")

for i, tile in enumerate(tiles):
    r = i // cols
    c = i % cols
    x = c * tw
    y = r * th
    collage.paste(tile, (x, y))

collage.save(args.output, quality=90)
