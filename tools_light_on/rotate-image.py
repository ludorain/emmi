#! /usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import skimage

import sys
import os

current_dir = os.path.dirname(__file__)
libs_path = os.path.abspath(os.path.join(current_dir, '..', 'libs'))
sys.path.insert(0, libs_path)

import emmi


def parse_arguments():
    parser = argparse.ArgumentParser(description='Rotate an EMMI image')
    parser.add_argument('--input', type=str, required=True, help='Input TIF filename')
    parser.add_argument('--angle', type=float, required=True, help='Rotation angle in degrees')
    parser.add_argument('--output', type=str, required=True, help='Output TIF filename')
    parser.add_argument('--display', action='store_true', help='Display original and rotated image')
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
    ax2.set_title('rotated')
    ax2.axis('on')
    im2 = ax2.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(im2, ax=ax2)

    plt.tight_layout()
    plt.show()


def rotated_valid_mask(image, angle):
    mask = np.ones(image.shape, dtype=float)
    mask = skimage.transform.rotate(mask, angle, resize=False, order=0)
    return mask > 0


def count_inserted_pixels(mask):
    return np.count_nonzero(~mask)


def get_clip_margins(mask):
    heights = np.zeros(mask.shape[1], dtype=int)
    best_area = 0
    best_rect = (0, mask.shape[0], 0, mask.shape[1])

    for row, values in enumerate(mask):
        heights = np.where(values, heights + 1, 0)
        stack = []

        for col in range(mask.shape[1] + 1):
            height = heights[col] if col < mask.shape[1] else 0

            while stack and heights[stack[-1]] > height:
                top_height = heights[stack.pop()]
                left = stack[-1] + 1 if stack else 0
                right = col
                area = top_height * (right - left)

                if area > best_area:
                    best_area = area
                    best_rect = (row - top_height + 1, row + 1, left, right)

            stack.append(col)

    top, bottom, left, right = best_rect
    return {
        'top': top,
        'bottom': mask.shape[0] - bottom,
        'left': left,
        'right': mask.shape[1] - right,
    }


if __name__ == "__main__":
    args = parse_arguments()

    print(' --- opening input image:', args.input)
    image = emmi.io.read_image(args.input)
    original = np.copy(image)

    image = emmi.rotation.rotate_image(image, args.angle)
    valid_mask = rotated_valid_mask(original, args.angle)
    inserted_pixels = count_inserted_pixels(valid_mask)
    inserted_fraction = inserted_pixels / original.size
    print(f' --- inserted black pixels: {inserted_pixels} ({inserted_fraction:.2%})')
    margins = get_clip_margins(valid_mask)
    print(' --- clip margins to remove black pixels: '
          f'top={margins["top"]}, bottom={margins["bottom"]}, '
          f'left={margins["left"]}, right={margins["right"]}')

    print(' --- saving output image:', args.output)
    emmi.io.save_image(args.output, image)

    if args.display:
        display(original, image, args.cmap, args.vrange)