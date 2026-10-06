#!/usr/bin/env python3

# Code to find luminous centers of a TIF image.
# Extended with:
#   --output_focus BASENAME     save focused crop as BASENAME.png and BASENAME.pdf
#   --focus_xy X Y              focus on explicit image coordinates, independently of source detection

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import tifffile

from photutils.background import Background2D, MedianBackground
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from photutils.segmentation import detect_sources
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.segmentation import SourceCatalog
from matplotlib.patches import Circle


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process an EMMI image')
    parser.add_argument('--input', type=str, required=True, help='Input TIF filename')
    parser.add_argument('--output_origin', type=str, required=False, help='Output PNG filename')
    parser.add_argument('--output', type=str, required=False, help='Output processed PNG filename')
    parser.add_argument('--display_original', action='store_true', help='Display image')
    parser.add_argument('--display_processed', action='store_true', help='Display processed image')
    parser.add_argument('--convolution', action='store_true', help='Apply convolution to the image')
    parser.add_argument('--coordinates_python', type=str, required=False, help='Output TXT filename for coordinates')
    parser.add_argument('--coordinates_root', type=str, required=False, help='Output TXT filename for coordinates')

    focus_group = parser.add_mutually_exclusive_group()
    focus_group.add_argument(
        '--focus',
        type=int,
        required=False,
        help='Show a region around the selected defect number (1-based index)'
    )
    focus_group.add_argument(
        '--focus_xy',
        type=float,
        nargs=2,
        metavar=('X', 'Y'),
        required=False,
        help='Show a region around explicit image coordinates X Y'
    )

    parser.add_argument(
        '--output_focus',
        type=str,
        required=False,
        help='Output basename for focused image; both PNG and PDF are saved'
    )
    parser.add_argument(
        '--focus_title',
        type=str,
        required=False,
        help='Custom title for the focused crop'
    )
    parser.add_argument('--circles', action='store_true', help='Display background-subtracted image with circles on detected sources')
    return parser.parse_args()


def make_focus_plot(processed, xc, yc, title, output_focus=None):
    """Crop a 100x100-pixel region around (xc, yc), draw the target circle and save/show it."""
    half_size = 50
    ny, nx = processed.shape

    if not (np.isfinite(xc) and np.isfinite(yc)):
        raise ValueError(f'focus coordinates are not finite: x={xc}, y={yc}')

    if xc < 0 or xc >= nx or yc < 0 or yc >= ny:
        raise ValueError(
            f'focus coordinates are outside the image: x={xc:.2f}, y={yc:.2f}, '
            f'image size={nx}x{ny}'
        )

    x_min = max(0, int(round(xc)) - half_size)
    x_max = min(nx, int(round(xc)) + half_size)
    y_min = max(0, int(round(yc)) - half_size)
    y_max = min(ny, int(round(yc)) + half_size)

    focus_image = processed[y_min:y_max, x_min:x_max]

    xc_local = xc - x_min
    yc_local = yc - y_min

    norm = ImageNormalize(focus_image, stretch=SqrtStretch())

    fig, ax = plt.subplots(figsize=(7.4, 7.6))
    ax.imshow(focus_image, origin='upper', norm=norm)

    circle = Circle(
        (xc_local, yc_local),
        radius=20,
        edgecolor='red',
        facecolor='none',
        linewidth=1.5
    )
    ax.add_patch(circle)

    ax.set_title(title, fontsize=14, pad=14, wrap=True)
    ax.set_xlim(0, focus_image.shape[1])
    ax.set_ylim(focus_image.shape[0], 0)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if output_focus:
        png_name = output_focus + '.png'
        pdf_name = output_focus + '.pdf'

        print(' --- saving focused image:', png_name)
        print(' --- saving focused image:', pdf_name)

        fig.savefig(png_name, dpi=300, bbox_inches='tight')
        fig.savefig(pdf_name, bbox_inches='tight')
    else:
        plt.show()

    plt.close(fig)


if __name__ == '__main__':
    args = parse_arguments()

    # Read the TIF image
    print(' --- opening input image:', args.input)
    image = tifffile.imread(args.input)

    # Display configuration
    plt.figure(figsize=(10, 5))
    plt.imshow(image)
    plt.axis('off')

    if args.output_origin:
        plt.savefig(args.output_origin, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
        plt.close()
    if args.display_original:
        plt.show()

    # Background subtraction
    bkg_estimator = MedianBackground()
    bkg = Background2D(
        image,
        (50, 50),
        filter_size=(3, 3),
        bkg_estimator=bkg_estimator
    )

    # Keep a floating-point working copy so background subtraction is always safe.
    processed = image.astype(float, copy=True)
    processed -= bkg.background

    # Define the detection threshold
    threshold = 3.0 * bkg.background_rms

    # Convolve the data with a 2D Gaussian kernel with FWHM of 3 pixels
    kernel = make_2dgaussian_kernel(3.0, size=5)
    convolved_data = convolve(processed, kernel)

    if not args.convolution:
        convolved_data = processed

    # Detect sources
    segment_map = detect_sources(convolved_data, threshold, npixels=30)
    print(segment_map)

    # Build source catalogue when at least one source is detected.
    cat = None
    tbl = None
    if segment_map is not None:
        cat = SourceCatalog(processed, segment_map, convolved_data=processed)
        print(cat)

        tbl = cat.to_table()
        tbl['xcentroid'].info.format = '.2f'
        tbl['ycentroid'].info.format = '.2f'
        print(tbl)
    else:
        print(' --- no sources detected in this image')

    # Processed images display and saving
    if args.display_processed or args.output:
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))

        ax1.imshow(processed, origin='upper', norm=norm)
        ax1.set_title('Background-subtracted Data')

        if segment_map is not None:
            ax2.imshow(
                segment_map,
                origin='upper',
                cmap=segment_map.cmap,
                interpolation='nearest'
            )
            ax2.set_title('Segmentation Image')
        else:
            ax2.text(0.5, 0.5, 'No sources detected', ha='center', va='center')
            ax2.set_axis_off()

        plt.tight_layout()

        if args.output:
            plt.savefig(args.output, dpi=300, bbox_inches='tight')

        if args.display_processed:
            plt.show()

        plt.close(fig)

    # Focus on a specific detected defect, preserving the original behaviour.
    if args.focus is not None:
        if tbl is None:
            print('Error: no defects were detected, therefore --focus cannot be used.')
            sys.exit(1)

        focus_index = args.focus - 1

        if focus_index < 0 or focus_index >= len(tbl):
            print(
                f'Error: selected defect {args.focus} is out of range. '
                f'Detected defects: {len(tbl)}'
            )
            sys.exit(1)

        xc = float(tbl['xcentroid'][focus_index])
        yc = float(tbl['ycentroid'][focus_index])

        print(f' --- focus on defect {args.focus}: x={xc:.2f}, y={yc:.2f}')

        try:
            make_focus_plot(
                processed,
                xc,
                yc,
                title=(args.focus_title if args.focus_title else f'Focus on defect {args.focus}'),
                output_focus=args.output_focus
            )
        except ValueError as exc:
            print(f'Error: {exc}')
            sys.exit(1)

    # Focus directly on supplied coordinates. This mode does not require the
    # hotspot to be detected in the current image.
    if args.focus_xy is not None:
        xc, yc = args.focus_xy

        print(f' --- focus on fixed coordinates: x={xc:.2f}, y={yc:.2f}')

        try:
            make_focus_plot(
                processed,
                xc,
                yc,
                title=(args.focus_title if args.focus_title else f'Focus at x={xc:.2f}, y={yc:.2f}'),
                output_focus=args.output_focus
            )
        except ValueError as exc:
            print(f'Error: {exc}')
            sys.exit(1)

    # Save .txt file with the coordinates of the centers
    if args.coordinates_python:
        if tbl is None:
            print('WARNING: no detected sources; writing an empty coordinate file.')
        print(' --- saving coordinates to:', args.coordinates_python)
        with open(args.coordinates_python, 'w') as f:
            if tbl is not None:
                for row in tbl:
                    f.write(f"{row['xcentroid']:.2f}, {row['ycentroid']:.2f}\n")

    if args.coordinates_root:
        if tbl is None:
            print('WARNING: no detected sources; writing an empty coordinate file.')
        print(' --- saving coordinates to:', args.coordinates_root)
        with open(args.coordinates_root, 'w') as f:
            if tbl is not None:
                image_height = processed.shape[0]
                for row in tbl:
                    x = row['xcentroid']
                    y = image_height - row['ycentroid']
                    f.write(f'{x:.2f}, {y:.2f}\n')

    # Display background-subtracted image with circles on detected sources
    if args.circles:
        if tbl is None:
            print('WARNING: --circles requested, but no sources were detected.')
        else:
            norm = ImageNormalize(stretch=SqrtStretch())

            fig, ax = plt.subplots(figsize=(8, 8))
            ax.imshow(processed, origin='upper', norm=norm)
            ax.set_title('Background-subtracted Data with detected sources')

            for row in tbl:
                x = row['xcentroid']
                y = row['ycentroid']

                circle = Circle(
                    (x, y),
                    radius=20,
                    edgecolor='red',
                    facecolor='none',
                    linewidth=1.5
                )
                ax.add_patch(circle)

            ax.set_xlim(0, processed.shape[1])
            ax.set_ylim(processed.shape[0], 0)

            plt.tight_layout()
            plt.show()
            plt.close(fig)
