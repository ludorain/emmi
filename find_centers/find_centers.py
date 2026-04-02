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
    parser.add_argument('--output_origin', type=str, required=False, help='Output PNG filename')
    parser.add_argument('--output', type=str, required=False, help='Output processed PNG filename')
    parser.add_argument('--display_original', action='store_true', help='Display image')
    parser.add_argument('--display_processed', action='store_true', help='Display processed image')
    parser.add_argument('--convolution', action='store_true', help='Apply convolution to the image')
    parser.add_argument('--coordinates', type=str, required=False, help='Output TXT filename for coordinates')
    parser.add_argument('--focus', type=int, required=False, help='Show a region around the selected defect number (1-based index)')
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

    if args.output_origin:
        plt.savefig(args.output_origin, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
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
        
    #Sources properties
    cat = SourceCatalog(processed, segment_map, convolved_data=processed)
    print(cat)

    tbl = cat.to_table()
    tbl['xcentroid'].info.format = '.2f'  # optional format
    tbl['ycentroid'].info.format = '.2f'
    #tbl['kron_flux'].info.format = '.2f'
    print(tbl)


    #Processed images display and saving
    if args.display_processed or args.output:
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))

        ax1.imshow(processed, origin='upper', norm=norm)
        ax1.set_title('Background-subtracted Data')

        ax2.imshow(segment_map, origin='upper', cmap=segment_map.cmap,
                interpolation='nearest')
        ax2.set_title('Segmentation Image')

        plt.tight_layout()

        if args.output:
            plt.savefig(args.output, dpi=300, bbox_inches='tight')

        if args.display_processed:
            plt.show()

        plt.close(fig)


    #Focus on a specific defect if requested
    if args.focus is not None:
        focus_index = args.focus - 1  # convert from 1-based to 0-based index

        if focus_index < 0 or focus_index >= len(tbl):
            print(f"Error: selected defect {args.focus} is out of range. "
                f"Detected defects: {len(tbl)}")
            exit()

        # Centroid coordinates of the selected defect
        xc = tbl['xcentroid'][focus_index]
        yc = tbl['ycentroid'][focus_index]

        print(f" --- focus on defect {args.focus}: x={xc:.2f}, y={yc:.2f}")

        # Half-size of the crop: 50 pixels -> total region about 100x100
        half_size = 50

        # Image shape
        ny, nx = processed.shape

        # Crop boundaries, clipped to image edges
        x_min = max(0, int(round(xc)) - half_size)
        x_max = min(nx, int(round(xc)) + half_size)
        y_min = max(0, int(round(yc)) - half_size)
        y_max = min(ny, int(round(yc)) + half_size)

        # Crop image
        focus_image = processed[y_min:y_max, x_min:x_max]

        # Coordinates of the centroid inside the cropped image
        xc_local = xc - x_min
        yc_local = yc - y_min

        # Display focused image
        norm = ImageNormalize(focus_image, stretch=SqrtStretch())

        plt.figure(figsize=(6, 6))
        plt.imshow(focus_image, origin='upper', norm=norm)
        plt.plot(xc_local, yc_local, 'rx', markersize=12, markeredgewidth=2)
        plt.title(f'Focus on defect {args.focus}')
        plt.xlim(0, focus_image.shape[1])
        plt.ylim(focus_image.shape[0], 0)
        plt.tight_layout()
        plt.show()

    #Save .txt file with the coordinates of the centers

    if args.coordinates:
        print(' --- saving coordinates to:', args.coordinates)
        with open(args.coordinates, 'w') as f:
            for row in tbl:
                f.write(f"{row['xcentroid']:.2f}, {row['ycentroid']:.2f}\n")
                

