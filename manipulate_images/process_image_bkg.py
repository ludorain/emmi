#! /usr/bin/env python
import sys
import os

#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import matplotlib.pyplot as plt
import numpy as np
import argparse
import skimage

current_dir = os.path.dirname(__file__)
libs_path = os.path.abspath(os.path.join(current_dir, '..', 'libs'))
sys.path.insert(0, libs_path)

import emmi

from photutils.background import Background2D, MedianBackground


def remove_background(image):
    bkg_estimator = MedianBackground()
    bkg = Background2D(
        image,
        (50, 50),
        filter_size=(3, 3),
        bkg_estimator=bkg_estimator
    )
    processed = np.copy(image)
    processed -= bkg.background
    return processed


process_dictionary = {
    'remove_column_bias'  : lambda image: emmi.process.remove_column_bias(image),
    'remove_hot_pixels'   : lambda image: emmi.process.remove_hot_pixels(image),
    'remove_cold_pixels'  : lambda image: emmi.process.remove_cold_pixels(image),
    'fix_underflows'      : lambda image: emmi.process.fix_underflows(image),
    'rescale_intensity'   : lambda image: emmi.process.rescale_intensity(image),
    'restore_intensity'   : lambda image: emmi.process.restore_intensity(image),
    'contrast_stretching' : lambda image: emmi.process.contrast_stretching(image),
    'equalize_adapthist'  : lambda image: emmi.process.equalize_adapthist(image),
    'equalize_hist'       : lambda image: emmi.process.equalize_hist(image),
    'threshold_local'     : lambda image: emmi.process.threshold_local(image),
    'threshold_otsu'      : lambda image: emmi.process.threshold_otsu(image),
    'denoise_nl_means'    : lambda image: emmi.process.denoise_nl_means(image),
    'sobel'               : lambda image: emmi.process.sobel(image),
    'canny'               : lambda image: emmi.process.canny(image),
    'remove_background'   : lambda image: remove_background(image),
}

process_available = list(process_dictionary.keys())


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process an EMMI image')
    parser.add_argument('--input', type=str, required=True, help='Input TIF filename')
    parser.add_argument('--output', type=str, required=False, help='Output TIF filename')
    parser.add_argument('--print', type=str, required=False, help='Print PNG filename')
    parser.add_argument('--process', type=str, nargs='+', required=False, choices=process_available, help='Process image')
    parser.add_argument('--display', action='store_true', help='Display original and processed image')
    parser.add_argument('--cmap', type=str, required=False, default=None, help='Display with requested cmap')
    parser.add_argument('--vrange', type=float, nargs=2, required=False, default=(None,None), help='Display with requested vmin/vmax')
    return parser.parse_args()


def display(original, image, cmap, vrange):
    fig = plt.figure(figsize=(20, 10))

    vmin, vmax = vrange

    ax1 = fig.add_subplot(2, 2, 1)
    ax1.set_title('original')
    ax1.axis('on')
    im1 = ax1.imshow(original, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(im1, ax=ax1)

    ax2 = fig.add_subplot(2, 2, 2, sharex=ax1, sharey=ax1)
    ax2.set_title('processed')
    ax2.axis('on')
    im2 = ax2.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(im2, ax=ax2)

    ax3 = fig.add_subplot(2, 2, 3)
    x, y = skimage.exposure.histogram(original)
    x, y = skimage.exposure.cumulative_distribution(original)
    ax3.plot(y, x)
    ax3.set_title('original cumulative distribution')

    ax4 = fig.add_subplot(2, 2, 4)
    x, y = skimage.exposure.histogram(image)
    x, y = skimage.exposure.cumulative_distribution(image)
    ax4.plot(y, x)
    ax4.set_title('processed cumulative distribution')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    args = parse_arguments()

    print(' --- opening input image:', args.input)
    image = emmi.io.read_image(args.input)
    original = np.copy(image)

    if args.process is not None:
        print(' --- process image')
        for process in args.process:
            print(f'     -> applying {process}')
            image = process_dictionary[process](image)

    if args.output is not None:
        print(' --- saving output image:', args.output)
        emmi.io.save_image(args.output, image)

    if args.print is not None:
        print(' --- print PNG image:', args.print)
        vmin, vmax = args.vrange
        emmi.io.print_png(args.print, image, args.cmap, vmin, vmax)

    if args.display:
        display(original, image, args.cmap, args.vrange)