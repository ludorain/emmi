#Code to find luminous centers of TIF image
#to compile, first: conda activate astropy
#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import argparse
import skimage
import tifffile


from photutils.background import Background2D, MedianBackground
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from photutils.segmentation import detect_sources
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.segmentation import SourceCatalog


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process an EMMI image')
    parser.add_argument('--input', type=str, required=True, help='Input TIF filename')
    parser.add_argument('--output', type=str, required=False, help='Output PNG filename')
    parser.add_argument('--display_original', action='store_true', help='Display image')
    parser.add_argument('--display_processed', action='store_true', help='Display processed image')
    parser.add_argument('--convolution', action='store_true', help='Apply convolution to the image')
    parser.add_argument('--coordinates', type=str, required=False, help='Output TXT filename for coordinates')
    return parser.parse_args()    

if __name__ == "__main__":
    
    args = parse_arguments()

    #Read the TIF image
    print(' --- opening input image:', args.input)
    image = tifffile.imread(args.input)

    
    #Display configuration 
    plt.figure(figsize=(10, 5))
    plt.imshow(image)
    plt.axis('off')

    if args.output:
        plt.savefig(args.output, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
        plt.close()
    if args.display_original:
        plt.show()

    #Background subtraction --magari metterlo opzionale
    
    bkg_estimator = MedianBackground()
    bkg = Background2D(image, (50, 50), filter_size=(3, 3),
                    bkg_estimator=bkg_estimator)
    processed = image               
    processed -= bkg.background  # subtract the background


    #Define the detection threshold
    threshold = 2.5 * bkg.background_rms

    #Convolve the data with a 2D Gaussian kernel with FWHM of 3 pixels
    kernel = make_2dgaussian_kernel(3.0, size=5)  # FWHM = 3.0
    convolved_data = convolve(processed, kernel)

    if not args.convolution:
        convolved_data = processed

    #Detect sources (eventually with convolution)
    segment_map = detect_sources(convolved_data, threshold, npixels=25)
    print(segment_map)

    if args.display_processed:
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
        ax1.imshow(processed, origin='upper', norm=norm)
        ax1.set_title('Background-subtracted Data')
        ax2.imshow(segment_map, origin='upper', cmap=segment_map.cmap,
                interpolation='nearest')
        ax2.set_title('Segmentation Image')
        plt.tight_layout()
        plt.show()


    #Sources properties
    cat = SourceCatalog(processed, segment_map, convolved_data=processed)
    print(cat)

    tbl = cat.to_table()
    tbl['xcentroid'].info.format = '.2f'  # optional format
    tbl['ycentroid'].info.format = '.2f'
    #tbl['kron_flux'].info.format = '.2f'
    print(tbl)


    #Save .txt file with the coordinates of the centers

    if args.coordinates:
        print(' --- saving coordinates to:', args.coordinates)
        with open(args.coordinates, 'w') as f:
            for row in tbl:
                f.write(f"{row['xcentroid']:.2f}, {row['ycentroid']:.2f}\n")
                

