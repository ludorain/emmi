
# python find_defects_isolated.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/2processed/run=20260513-032137_x=0_y=0_z=0_T=20_v=58.30_data=diff_processed.tif" --circles

#! /usr/bin/env python

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import tifffile

from astropy.convolution import convolve
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from matplotlib.patches import Circle
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import SourceCatalog, SourceFinder, make_2dgaussian_kernel


# ============================================================
# USER CONFIGURATION
# ============================================================
# Radius used later to integrate the luminosity around each
# isolated source. It is NOT estimated from the segmented source.
integration_radius = 8.0  # pixels

# Minimum required distance between the centroid of an isolated
# source and the centroid of every other detected source.
isolation_radius = 60.0  # pixels


# Fixed display scale used by --focus and --focus_area.
vmin_global = 0.0
vmax_global = 20.0


# ============================================================
# ARGUMENTS
# ============================================================
def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Detect and deblend EMMI sources, then select only sources '
            'isolated from every other detected source.'
        )
    )

    parser.add_argument('--input', type=str, required=True,
                        help='Input TIF filename')
    parser.add_argument('--output_origin', type=str, required=False,
                        help='Output PNG filename for the original image')
    parser.add_argument('--output', type=str, required=False,
                        help='Output PNG filename for the processed image')

    parser.add_argument('--display_original', action='store_true',
                        help='Display original image')
    parser.add_argument('--display_processed', action='store_true',
                        help='Display processed image and segmentation map')

    parser.add_argument('--coordinates_python', type=str, required=False,
                        help='Output TXT filename for isolated-source coordinates')
    parser.add_argument('--coordinates_root', type=str, required=False,
                        help='Output TXT filename for isolated-source ROOT coordinates')

    parser.add_argument('--focus', type=int, required=False,
                        help=(
                            'Show a region around an isolated defect_id written '
                            'in --coordinates_python'
                        ))

    parser.add_argument('--focus_area', nargs=3, type=float,
                        metavar=('X', 'Y', 'RADIUS'),
                        help='Focus on a manually selected area: x y radius')

    parser.add_argument('--circles', action='store_true',
                        help=(
                            'Display all detected centroids and draw integration '
                            'circles around isolated sources'
                        ))

    # SourceFinder / deblending parameters
    parser.add_argument('--npixels', type=int, default=30,
                        help='Minimum number of connected pixels for source detection')
    parser.add_argument('--nlevels', type=int, default=32,
                        help='Number of multi-thresholding levels for deblending')
    parser.add_argument('--contrast', type=float, default=0.001,
                        help='Deblending contrast parameter')

    return parser.parse_args()


# ============================================================
# HELPER FUNCTIONS
# ============================================================
def validate_radii():
    """Validate the fixed integration and isolation radii."""
    if integration_radius <= 0:
        raise ValueError(
            f'Invalid integration_radius={integration_radius}. '
            'integration_radius must be greater than zero.'
        )

    if isolation_radius <= 2.0 * integration_radius:
        raise ValueError(
            'Invalid radius configuration: isolation_radius must be greater '
            'than 2 * integration_radius. '
            f'Current values: isolation_radius={isolation_radius}, '
            f'integration_radius={integration_radius}, '
            f'2*integration_radius={2.0 * integration_radius}.'
        )


def get_centroid_column_names(tbl):
    """Return centroid column names for different Photutils versions."""
    if 'x_centroid' in tbl.colnames and 'y_centroid' in tbl.colnames:
        return 'x_centroid', 'y_centroid'

    if 'xcentroid' in tbl.colnames and 'ycentroid' in tbl.colnames:
        return 'xcentroid', 'ycentroid'

    raise RuntimeError('Centroid columns not found in SourceCatalog table.')


def to_float(value):
    """Convert an Astropy Quantity or table value to a plain float."""
    if hasattr(value, 'value'):
        value = value.value
    return float(value)


def select_isolated_sources(tbl, x_col, y_col):
    """
    Select sources whose centroid is farther than isolation_radius from
    every other detected source.

    Returns
    -------
    isolated_tbl : astropy.table.Table
        Sub-table containing only isolated sources.
    nearest_distances : numpy.ndarray
        Nearest-neighbor distance for every detected source.
    isolated_mask : numpy.ndarray
        Boolean mask relative to the complete source table.
    """
    x = np.array([to_float(value) for value in tbl[x_col]], dtype=float)
    y = np.array([to_float(value) for value in tbl[y_col]], dtype=float)

    coordinates = np.column_stack((x, y))
    number_of_sources = len(coordinates)

    if number_of_sources == 1:
        nearest_distances = np.array([np.inf])
        isolated_mask = np.array([True])
    else:
        differences = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
        distance_matrix = np.sqrt(np.sum(differences**2, axis=2))

        # Exclude each source's zero distance from itself.
        np.fill_diagonal(distance_matrix, np.inf)

        nearest_distances = np.min(distance_matrix, axis=1)
        isolated_mask = nearest_distances > isolation_radius

    isolated_tbl = tbl[isolated_mask].copy()

    # defect_id refers only to the isolated-source list written by
    # --coordinates_python. It is 1-based and consecutive.
    isolated_tbl['defect_id'] = np.arange(1, len(isolated_tbl) + 1)
    isolated_tbl['nearest_source_distance_pix'] = nearest_distances[isolated_mask]

    return isolated_tbl, nearest_distances, isolated_mask


def display_focus(processed, xc, yc, radius, title):
    """Display a cropped image around a selected point."""
    half_size = max(20, int(np.ceil(2.0 * radius)))

    ny, nx = processed.shape

    x_min = max(0, int(round(xc)) - half_size)
    x_max = min(nx, int(round(xc)) + half_size)
    y_min = max(0, int(round(yc)) - half_size)
    y_max = min(ny, int(round(yc)) + half_size)

    focus_image = processed[y_min:y_max, x_min:x_max]

    if focus_image.size == 0:
        raise ValueError(
            f'The requested focus area around x={xc}, y={yc} is outside the image.'
        )

    xc_local = xc - x_min
    yc_local = yc - y_min

    norm = ImageNormalize(
        vmin=vmin_global,
        vmax=vmax_global,
        stretch=SqrtStretch()
    )

    fig_focus, ax_focus = plt.subplots(figsize=(6, 6))
    ax_focus.imshow(focus_image, origin='upper', norm=norm)

    circle = Circle(
        (xc_local, yc_local),
        radius=radius,
        edgecolor='red',
        facecolor='none',
        linewidth=1.5
    )
    ax_focus.add_patch(circle)

    ax_focus.set_title(title)
    ax_focus.set_xlim(0, focus_image.shape[1])
    ax_focus.set_ylim(focus_image.shape[0], 0)

    plt.tight_layout()
    plt.show()
    plt.close(fig_focus)


# ============================================================
# MAIN PROGRAM
# ============================================================
def main():
    args = parse_arguments()

    try:
        validate_radii()
    except ValueError as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1

    fixed_area = np.pi * integration_radius**2

    print(f' --- integration_radius = {integration_radius:.2f} pixels')
    print(f' --- integration area   = {fixed_area:.2f} pixels^2')
    print(f' --- isolation_radius   = {isolation_radius:.2f} pixels')

    # ------------------------------------------------------------
    # Read the TIF image
    # ------------------------------------------------------------
    print(' --- opening input image:', args.input)
    image = tifffile.imread(args.input).astype(float)

    if image.ndim != 2:
        print(
            f'ERROR: expected a 2D image, but input shape is {image.shape}.',
            file=sys.stderr
        )
        return 1

    # ------------------------------------------------------------
    # Display / save original image
    # ------------------------------------------------------------
    if args.display_original or args.output_origin:
        fig0, ax0 = plt.subplots(figsize=(10, 5))
        ax0.imshow(image, origin='upper')
        ax0.axis('off')

        if args.output_origin:
            fig0.savefig(
                args.output_origin,
                format='png',
                dpi=300,
                bbox_inches='tight',
                pad_inches=0
            )

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
    # Convolution: always performed and always used for detection
    # ------------------------------------------------------------
    kernel = make_2dgaussian_kernel(3.0, size=5)  # FWHM = 3.0 pixels
    convolved_data = convolve(processed, kernel)
    detection_image = convolved_data

    # ------------------------------------------------------------
    # Detection and deblending of ALL sources
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
        print(' --- No sources detected.')
        return 0

    # ------------------------------------------------------------
    # Centroids only: no segmented area or equivalent radius is used
    # ------------------------------------------------------------
    cat = SourceCatalog(
        processed,
        segment_map,
        convolved_data=detection_image
    )

    tbl = cat.to_table()
    x_col, y_col = get_centroid_column_names(tbl)

    tbl[x_col].info.format = '.2f'
    tbl[y_col].info.format = '.2f'

    # ID in the complete detected/deblended source list.
    tbl['all_source_id'] = np.arange(1, len(tbl) + 1)

    # ------------------------------------------------------------
    # Select the isolated-source sub-category
    # ------------------------------------------------------------
    isolated_tbl, nearest_distances, isolated_mask = select_isolated_sources(
        tbl,
        x_col,
        y_col
    )

    number_all = len(tbl)
    number_isolated = len(isolated_tbl)
    number_non_isolated = number_all - number_isolated

    print(f' --- detected and deblended sources: {number_all}')
    print(f' --- isolated sources:              {number_isolated}')
    print(f' --- non-isolated sources:          {number_non_isolated}')

    if number_isolated > 0:
        print(' --- isolated-source centroids:')
        for row in isolated_tbl:
            print(
                f"     defect_id={int(row['defect_id'])}, "
                f"x={to_float(row[x_col]):.2f}, "
                f"y={to_float(row[y_col]):.2f}, "
                f"nearest distance="
                f"{to_float(row['nearest_source_distance_pix']):.2f} pixels"
            )

    # ------------------------------------------------------------
    # Processed image and deblended segmentation map
    # ------------------------------------------------------------
    if args.display_processed or args.output:
        norm = ImageNormalize(stretch=SqrtStretch())

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))

        ax1.imshow(processed, origin='upper', norm=norm)
        ax1.set_title('Background-subtracted Data')
        ax1.axis('off')

        ax2.imshow(
            segment_map.data,
            origin='upper',
            cmap=segment_map.cmap,
            interpolation='nearest'
        )
        ax2.set_title('Deblended Segmentation Image — All Sources')
        ax2.axis('off')

        plt.tight_layout()

        if args.output:
            fig.savefig(args.output, dpi=300, bbox_inches='tight')

        if args.display_processed:
            plt.show()

        plt.close(fig)

    # ------------------------------------------------------------
    # Focus on an isolated source
    # ------------------------------------------------------------
    if args.focus is not None:
        matches = np.where(
            np.array(isolated_tbl['defect_id'], dtype=int) == args.focus
        )[0]

        if len(matches) == 0:
            if number_isolated == 0:
                print('ERROR: no isolated sources are available for --focus.')
            else:
                print(
                    f'ERROR: selected isolated defect_id {args.focus} is out of '
                    f'range. Available values: 1 ... {number_isolated}'
                )
            return 1

        focus_index = int(matches[0])
        xc = to_float(isolated_tbl[x_col][focus_index])
        yc = to_float(isolated_tbl[y_col][focus_index])
        nearest_distance = to_float(
            isolated_tbl['nearest_source_distance_pix'][focus_index]
        )

        print(
            f' --- focus on isolated defect {args.focus}: '
            f'x={xc:.2f}, y={yc:.2f}, '
            f'integration radius={integration_radius:.2f} pixels, '
            f'nearest-source distance={nearest_distance:.2f} pixels'
        )

        display_focus(
            processed,
            xc,
            yc,
            integration_radius,
            f'Isolated defect {args.focus}'
        )

    # ------------------------------------------------------------
    # Focus on a manually selected area
    # ------------------------------------------------------------
    if args.focus_area is not None:
        xc, yc, manual_radius = args.focus_area

        if manual_radius <= 0:
            print('ERROR: the --focus_area radius must be greater than zero.')
            return 1

        print(
            f' --- focus on manual area: '
            f'x={xc:.2f}, y={yc:.2f}, radius={manual_radius:.2f} pixels'
        )

        display_focus(
            processed,
            xc,
            yc,
            manual_radius,
            f'Manual focus: x={xc:.1f}, y={yc:.1f}, r={manual_radius:.1f} px'
        )

    # ------------------------------------------------------------
    # Save isolated coordinates using the FIXED integration radius
    # ------------------------------------------------------------
    # coordinates_python columns:
    # defect_id, x, y, integration_area_pix2, integration_radius_pix
    #
    # coordinates_root columns:
    # x, y_root, integration_area_pix2, integration_radius_pix
    #
    # The area is pi * integration_radius^2. It is not obtained from
    # the segmented source footprint.
    if args.coordinates_python:
        print(' --- saving isolated coordinates to:', args.coordinates_python)

        with open(args.coordinates_python, 'w', encoding='utf-8') as output_file:
            for row in isolated_tbl:
                defect_id = int(row['defect_id'])
                x = to_float(row[x_col])
                y = to_float(row[y_col])

                output_file.write(
                    f'{defect_id}, {x:.2f}, {y:.2f}, '
                    f'{fixed_area:.2f}, {integration_radius:.2f}\n'
                )

    if args.coordinates_root:
        print(' --- saving isolated ROOT coordinates to:', args.coordinates_root)

        image_height = float(processed.shape[0])

        with open(args.coordinates_root, 'w', encoding='utf-8') as output_file:
            for row in isolated_tbl:
                x = to_float(row[x_col])
                y_root = image_height - to_float(row[y_col])

                output_file.write(
                    f'{x:.2f}, {y_root:.2f}, '
                    f'{fixed_area:.2f}, {integration_radius:.2f}\n'
                )

    # ------------------------------------------------------------
    # Display all centroids and isolated integration circles
    # ------------------------------------------------------------
    if args.circles:
        norm = ImageNormalize(stretch=SqrtStretch())

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.imshow(processed, origin='upper', norm=norm)
        ax.set_title(
            'Detected Sources and Isolated Integration Regions\n'
            f'integration radius = {integration_radius:.1f} px, '
            f'isolation radius = {isolation_radius:.1f} px'
        )

        # Mark every detected/deblended source centroid.
        for row in tbl:
            x = to_float(row[x_col])
            y = to_float(row[y_col])
            ax.plot(x, y, marker='x', markersize=4, color='yellow')

        # Draw fixed-radius integration circles only around isolated sources.
        for row in isolated_tbl:
            defect_id = int(row['defect_id'])
            x = to_float(row[x_col])
            y = to_float(row[y_col])

            circle = Circle(
                (x, y),
                radius=integration_radius,
                edgecolor='red',
                facecolor='none',
                linewidth=1.5
            )

            ax.add_patch(circle)
            ax.text(x, y, str(defect_id), color='red', fontsize=8)

        ax.set_xlim(0, processed.shape[1])
        ax.set_ylim(processed.shape[0], 0)

        plt.tight_layout()
        plt.show()
        plt.close(fig)

    return 0


if __name__ == '__main__':
    sys.exit(main())