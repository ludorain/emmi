#! /usr/bin/env python

import argparse
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import emmi

def display(image1, image2, stitched):

    fig = plt.figure(layout="constrained")
    gs = GridSpec(2, 2, figure=fig)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.imshow(image1, cmap='gray')
    ax1.set_title('image 1')
    ax1.axis('off')
    
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.imshow(image2, cmap='gray')
    ax2.set_title('image 2')
    ax2.axis('off')
    
    ax3 = fig.add_subplot(gs[1, :])
    ax3.imshow(stitched, cmap='gray')
    ax3.set_title('stitched')
    ax3.axis('off')
    
    plt.tight_layout()
    plt.show()    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="D")
    parser.add_argument("--input", type=str, nargs=2, required=True, help="Input TIF filenames to stitch")
    parser.add_argument("--output", required=False, help="Output TIF file name")
    parser.add_argument('--shift', type=float, nargs=2, default=None, help="Image offset (y, x)")
    parser.add_argument('--seed', type=float, nargs=2, default=None, help="Image offset seed (y, x)")
    parser.add_argument("--display", action="store_true", help="Display process")
    args = parser.parse_args()

    ### read images
    image1 = emmi.read_image(args.input[0])
    image2 = emmi.read_image(args.input[1])

    image1_mask = None
    image2_mask = None
    
    if False and args.seed is not None:
        image1_mask = np.zeros(image1.shape, dtype=bool)
        image2_mask = np.zeros(image2.shape, dtype=bool)
        if args.seed[1] > 0:
            image1_mask[ : , int(args.seed[1]):] = True
            image2_mask[ : , :-int(args.seed[1])] = True
        elif args.seed[1] < 0:
            image2_mask[ : , -int(args.seed[1]):] = True
            image1_mask[ : , :int(args.seed[1])] = True
        if args.seed[0] > 0:
            image1_mask[int(args.seed[0]): , :] = True
            image2_mask[:-int(args.seed[0]) , :] = True
        elif args.seed[0] < 0:
            image2_mask[-int(args.seed[0]): , :] = True
            image1_mask[:int(args.seed[0]) , :] = True
    

    if True and args.seed is not None:
        if args.seed[1] > 0:
            image1 = image1[ : , int(args.seed[1]):]
            image2 = image2[ : , :-int(args.seed[1])]
        elif args.seed[1] < 0:
            image2 = image2[ : , -int(args.seed[1]):]
            image1 = image1[ : , :int(args.seed[1])]

        if args.seed[0] > 0:
            image1 = image1[int(args.seed[0]): , :]
            image2 = image2[:-int(args.seed[0]) , :]
        elif args.seed[0] < 0:
            image2 = image2[-int(args.seed[0]): , :]
            image1 = image1[:int(args.seed[0]) , :]

    ### calculate shift
    if args.shift is None:
#        shift = emmi.phase_correlate(image1, image2)
        shift = emmi.phase_correlate(image1, image2, reference_mask=image1_mask, moving_mask=image2_mask)
    else:
        shift = args.shift

    ### stitch
    stitched = emmi.stitch_images(image1, image2, shift)
    
    ### write output
    if args.output is not None:
        emmi.save_image(args.output, stitched)

    ### display process
    if args.display:
        display(image1, image2, stitched)
