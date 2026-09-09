#! /usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np

import sys
import os

current_dir = os.path.dirname(__file__)
libs_path = os.path.abspath(os.path.join(current_dir, '..', 'libs'))
sys.path.insert(0, libs_path)

import emmi


def parse_arguments():
    parser = argparse.ArgumentParser(description='Shift an EMMI image')
    parser.add_argument('--input', type=str, required=True, help='Input TIF filename')
    parser.add_argument('--shift', type=float, nargs=2, required=True, help='Image shift as dy dx')
    parser.add_argument('--output', type=str, required=True, help='Output TIF filename')
    parser.add_argument('--display', action='store_true', help='Display original and shifted image')
    parser.add_argument('--cmap', type=str, required=False, default=None, help='Display with requested cmap')
    parser.add_argument('--vrange', type=float, nargs=2, required=False, default=(None,None), help='Display with requested vmin/vmax')
    return parser.parse_args()


def display(original, image, cmap, vrange):
    fig = plt.figure(figsize=(20, 10))

    vmin, vmax = vrange

    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_title('original')
    ax1.axis('on')
    im1 = ax1.imshow(original, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(im1, ax=ax1)

    ax2 = fig.add_subplot(1, 2, 2, sharex=ax1, sharey=ax1)
    ax2.set_title('shifted')
    ax2.axis('on')
    im2 = ax2.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(im2, ax=ax2)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    args = parse_arguments()

    print(' --- opening input image:', args.input)
    image = emmi.io.read_image(args.input)
    original = np.copy(image)

    image = emmi.stitching.shift_image(image, args.shift)

    print(' --- saving output image:', args.output)
    emmi.io.save_image(args.output, image)

    if args.display:
        display(original, image, args.cmap, args.vrange)