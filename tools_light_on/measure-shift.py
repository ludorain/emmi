#! /usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import emmi

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process an EMMI image')
    parser.add_argument('--input', type=str, nargs=2, required=True, help='Input TIF filenames')
    parser.add_argument('--shift', type=float, nargs=2, default=None, help="Shift between images")
    parser.add_argument('--display', action='store_true', help='Display original and processed image')
    parser.add_argument('--cmap', type=str, required=False, default=None, help='Display with requested cmap')
    parser.add_argument('--vrange', type=float, nargs=2, required=False, default=(None,None), help='Display with requested vmin/vmax')
    return parser.parse_args()    

def display(image1, image2, cmap, vrange):
    fig = plt.figure(figsize=(10, 5))
    vmin, vmax = vrange
    
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_title('image1')
    ax1.axis('on')
    im1 = ax1.imshow(image1, cmap=cmap, vmin=vmin, vmax=vmax)

    ax2 = fig.add_subplot(1, 2, 2, sharex=ax1, sharey=ax1)
    ax2.set_title('image2')
    ax2.axis('on')
    im2 = ax2.imshow(image2, cmap=cmap, vmin=vmin, vmax=vmax)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    args = parse_arguments()
    input1, input2 = args.input
    image1 = emmi.io.read_image(input1)
    image2 = emmi.io.read_image(input2)

    shift = args.shift
    if shift is None:
        shift = emmi.stitching.phase_correlate(image1, image2)
    image2 = emmi.stitching.shift_image(image2, shift)

    if args.display:
        display(image1, image2, args.cmap, args.vrange)
