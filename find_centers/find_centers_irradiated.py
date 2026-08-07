# Code to find luminous centers of TIF image
# This version main characteristics are:
# 1) The radius changes according to the dimension of the spot


# python find_centers_irradiated.py --input "/Users/ludovicarainero/emmi/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/1originals/run=20260513-032137_x=0_y=0_z=0_T=20_v=58.30_data=diff.tif" --circles
# python find_centers_irradiated.py --input "/Users/ludovicarainero/emmi/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/2processed/run=20260513-032137_x=0_y=0_z=0_T=20_v=53.30_data=diff_processed.tif" --convolution 
#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import argparse
import tifffile

from photutils.background import Background2D, MedianBackground
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from photutils.segmentation import SourceFinder
from photutils.segmentation import SourceCatalog
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from matplotlib.patches import Circle


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process an EMMI image')

    parser.add_argument('--input', type=str, required=True, help='Input TIF filename')
    parser.add_argument('--output_origin', type=str, required=False, help='Output PNG filename')
    parser.add_argument('--output', type=str, required=False, help='Output processed PNG filename')

    parser.add_argument('--display_original', action='store_true', help='Display original image')
    parser.add_argument('--display_processed', action='store_true', help='Display processed image')

    parser.add_argument('--convolution', action='store_true', help='Apply convolution to the image')

    parser.add_argument('--coordinates_python', type=str, required=False,
                        help='Output TXT filename for Python coordinates')
    parser.add_argument('--coordinates_root', type=str, required=False,
                        help='Output TXT filename for ROOT coordinates')

    parser.add_argument('--focus', type=int, required=False,
                        help='Show a region around the selected defect number, 1-based index')

    parser.add_argument('--circles', action='store_true',
                        help='Display background-subtracted image with circles on detected sources')

    # SourceFinder / deblending parameters
    parser.add_argument('--npixels', type=int, default=30,
                        help='Minimum number of connected pixels for source detection')
    parser.add_argument('--nlevels', type=int, default=32,
                        help='Number of multi-thresholding levels for deblending')
    parser.add_argument('--contrast', type=float, default=0.001,
                        help='Deblending contrast parameter')

    # Display parameters for circles
    parser.add_argument('--circle_scale', type=float, default=1.0,
                        help='Scale factor applied to equivalent radius when drawing circles')
    parser.add_argument('--circle_min_radius', type=float, default=3.0,
                        help='Minimum displayed circle radius in pixels')

    return parser.parse_args()


def get_centroid_column_names(tbl):
    """
    Photutils versions may use either:
        x_centroid, y_centroid
    or older names:
        xcentroid, ycentroid
    This function makes the rest of the code robust.
    """
    if 'x_centroid' in tbl.colnames and 'y_centroid' in tbl.colnames:
        return 'x_centroid', 'y_centroid'

    if 'xcentroid' in tbl.colnames and 'ycentroid' in tbl.colnames:
        return 'xcentroid', 'ycentroid'

    raise RuntimeError("Centroid columns not found in SourceCatalog table.")


def to_float(value):
    """
    Convert Astropy Quantity or table value to a plain float.
    """
    if hasattr(value, 'value'):
        value = value.value
    return float(value)


if __name__ == "__main__":

    args = parse_arguments()

    # ------------------------------------------------------------
    # Read the TIF image
    # ------------------------------------------------------------
    print(' --- opening input image:', args.input)

    # Convert to float to avoid problems when subtracting the background
    # from integer or unsigned images.
    image = tifffile.imread(args.input).astype(float)

    # ------------------------------------------------------------
    # Display / save original image
    # ------------------------------------------------------------
    if args.display_original or args.output_origin:
        fig0, ax0 = plt.subplots(figsize=(10, 5))
        ax0.imshow(image, origin='upper')
        ax0.axis('off')

        if args.output_origin:
            fig0.savefig(args.output_origin, format='png', dpi=300,
                         bbox_inches='tight', pad_inches=0)

        if args.display_original:
            plt.show()

        plt.close(fig0)

    # ------------------------------------------------------------
    # Background subtraction
    # ------------------------------------------------------------
    bkg_estimator = MedianBackground()

    bkg = Background2D(
        image,
        box_size=(50, 50),
        filter_size=(3, 3),
        bkg_estimator=bkg_estimator
    )

    processed = image - bkg.background

    # ------------------------------------------------------------
    # Detection threshold
    # ------------------------------------------------------------
    threshold = 4.0 * bkg.background_rms

    # ------------------------------------------------------------
    # Convolution
    # ------------------------------------------------------------
    kernel = make_2dgaussian_kernel(3.0, size=5)  # FWHM = 3.0 px
    convolved_data = convolve(processed, kernel)

    if args.convolution:
        detection_image = convolved_data
    else:
        detection_image = processed

    # ------------------------------------------------------------
    # Source detection + deblending using SourceFinder
    # ------------------------------------------------------------
    finder = SourceFinder(
        npixels=args.npixels,
        deblend=True,
        nlevels=args.nlevels,
        contrast=args.contrast,
        progress_bar=False
    )

    segment_map = finder(detection_image, threshold)

    if segment_map is None:
        print(" --- No sources detected.")
        exit()

    print(segment_map)

    # ------------------------------------------------------------
    # Source properties
    # ------------------------------------------------------------
    cat = SourceCatalog(
        processed,
        segment_map,
        convolved_data=detection_image
    )

    print(cat)

    tbl = cat.to_table()

    x_col, y_col = get_centroid_column_names(tbl)

    tbl[x_col].info.format = '.2f'
    tbl[y_col].info.format = '.2f'

    # ------------------------------------------------------------
    # Add source area and equivalent circular radius
    # ------------------------------------------------------------
    # segment_area is the area of the segmented source footprint.
    # equivalent_radius = sqrt(area/pi), assuming a circular source.
    segment_area = np.array([
        to_float(a) for a in cat.segment_area
    ])

    equivalent_radius = np.sqrt(segment_area / np.pi)

    tbl['segment_area_pix2'] = segment_area
    tbl['equivalent_radius_pix'] = equivalent_radius

    tbl['segment_area_pix2'].info.format = '.2f'
    tbl['equivalent_radius_pix'].info.format = '.2f'

    print(tbl)

    # ------------------------------------------------------------
    # Processed images display and saving
    # ------------------------------------------------------------
    if args.display_processed or args.output:
        norm = ImageNormalize(stretch=SqrtStretch())

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))

        ax1.imshow(processed, origin='upper', norm=norm)
        ax1.set_title('Background-subtracted Data')
        ax1.axis('off')

        ax2.imshow(segment_map.data, origin='upper',
                   cmap=segment_map.cmap, interpolation='nearest')
        ax2.set_title('Deblended Segmentation Image')
        ax2.axis('off')

        plt.tight_layout()

        if args.output:
            fig.savefig(args.output, dpi=300, bbox_inches='tight')

        if args.display_processed:
            plt.show()

        plt.close(fig)

    # ------------------------------------------------------------
    # Focus on a specific defect if requested
    # ------------------------------------------------------------
    if args.focus is not None:
        focus_index = args.focus - 1  # convert from 1-based to 0-based index

        if focus_index < 0 or focus_index >= len(tbl):
            print(f"Error: selected defect {args.focus} is out of range. "
                  f"Detected defects: {len(tbl)}")
            exit()

        # Centroid coordinates of the selected defect
        xc = to_float(tbl[x_col][focus_index])
        yc = to_float(tbl[y_col][focus_index])
        radius = to_float(tbl['equivalent_radius_pix'][focus_index])

        display_radius = max(
            args.circle_min_radius,
            args.circle_scale * radius
        )

        print(f" --- focus on defect {args.focus}: "
              f"x={xc:.2f}, y={yc:.2f}, "
              f"area={tbl['segment_area_pix2'][focus_index]:.2f} pix^2, "
              f"equivalent radius={radius:.2f} pix")

        # Half-size of the crop: 50 pixels -> total region about 100x100
        half_size = 20

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

        fig_focus, ax_focus = plt.subplots(figsize=(6, 6))
        ax_focus.imshow(focus_image, origin='upper', norm=norm)

        circle = Circle(
            (xc_local, yc_local),
            radius=display_radius,
            edgecolor='red',
            facecolor='none',
            linewidth=1.5
        )

        ax_focus.add_patch(circle)

        #ax_focus.set_title(f'Focus on defect {args.focus}')
        ax_focus.set_title('Focus on defect 2 - v_over = 4 V')
        ax_focus.set_xlim(0, focus_image.shape[1])
        ax_focus.set_ylim(focus_image.shape[0], 0)

        plt.tight_layout()
        plt.show()
        plt.close(fig_focus)

    # ------------------------------------------------------------
    # Save .txt file with coordinates + area + equivalent radius
    # ------------------------------------------------------------
    # Output columns:
    # x, y, segment_area_pix2, equivalent_radius_pix

    if args.coordinates_python:
        print(' --- saving coordinates to:', args.coordinates_python)

        with open(args.coordinates_python, 'w') as f:
            for row in tbl:
                x = to_float(row[x_col])
                y = to_float(row[y_col])
                area = to_float(row['segment_area_pix2'])
                radius = to_float(row['equivalent_radius_pix'])

                f.write(f"{x:.2f}, {y:.2f}, {area:.2f}, {radius:.2f}\n")

    if args.coordinates_root:
        print(' --- saving coordinates to:', args.coordinates_root)

        with open(args.coordinates_root, 'w') as f:
            for row in tbl:
                x = to_float(row[x_col])
                y = 1080.0 - to_float(row[y_col])
                area = to_float(row['segment_area_pix2'])
                radius = to_float(row['equivalent_radius_pix'])

                f.write(f"{x:.2f}, {y:.2f}, {area:.2f}, {radius:.2f}\n")

    # ------------------------------------------------------------
    # Display background-subtracted image with variable-size circles
    # ------------------------------------------------------------
    if args.circles:
        norm = ImageNormalize(stretch=SqrtStretch())

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.imshow(processed, origin='upper', norm=norm)
        ax.set_title('Background-subtracted Data with detected sources')

        for row in tbl:
            x = to_float(row[x_col])
            y = to_float(row[y_col])
            radius = to_float(row['equivalent_radius_pix'])

            display_radius = max(
                args.circle_min_radius,
                args.circle_scale * radius
            )

            circle = Circle(
                (x, y),
                radius=display_radius,
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