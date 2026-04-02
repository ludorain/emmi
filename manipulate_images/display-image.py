#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import argparse
import skimage

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process an EMMI image')
    parser.add_argument('--input', type=str, required=True, help='Input TIF filename')
    parser.add_argument('--output', type=str, required=False, help='Output PNG filename')
    parser.add_argument('--original', action='store_true', help='Display image')
    parser.add_argument('--display', action='store_true', help='Display image')
    parser.add_argument('--colorbar', action='store_true', help='Display image')
    parser.add_argument('--tightlayout', action='store_true', help='Display image')
    parser.add_argument('--cmap', type=str, required=False, help='')
    parser.add_argument('--vmin', type=float, required=False, help='')
    parser.add_argument('--vmax', type=float, required=False, help='')
    return parser.parse_args()    

def contrast_stretching(image):
    print(' --- contrast_stretching')
    pmin, pmax = np.percentile(image, (1, 99))
    image = skimage.exposure.rescale_intensity(image, in_range=(pmin, pmax))
    return image

if __name__ == "__main__":
    args = parse_arguments()

    print(' --- opening input image:', args.input)
    image = skimage.io.imread(args.input)
    if not args.original:
        image = contrast_stretching(image)
    
    plt.figure(figsize=(10, 5))
    plt.imshow(image, cmap=args.cmap, vmin=args.vmin, vmax=args.vmax)
    plt.axis('off')
    
    if args.colorbar:
        plt.colorbar()
    if args.tightlayout:
        plt.tight_layout()
    if args.output:
        plt.savefig(args.output, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
        plt.close()
    if args.display:
        plt.show()
