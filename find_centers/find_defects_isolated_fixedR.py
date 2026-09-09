#!/usr/bin/env python3

#For finding images lum vs overvoltage

#before annealing 
# v=1
# python find_defects_isolated_fixedR.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/2processed/run=20260513-032137_x=0_y=0_z=0_T=20_v=52.30_data=diff_processed.tif"   --focus_area 516 195 15

#v=2
# python find_defects_isolated_fixedR.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/2processed/run=20260513-032137_x=0_y=0_z=0_T=20_v=53.30_data=diff_processed.tif"   --focus_area 516 195 15

#v=3
# python find_defects_isolated_fixedR.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/2processed/run=20260513-032137_x=0_y=0_z=0_T=20_v=54.30_data=diff_processed.tif"   --focus_area 516 195 15

#v=4
# python find_defects_isolated_fixedR.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/2processed/run=20260513-032137_x=0_y=0_z=0_T=20_v=55.30_data=diff_processed.tif"   --focus_area 516 195 15


#v=5
# python find_defects_isolated_fixedR.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/2processed/run=20260513-032137_x=0_y=0_z=0_T=20_v=56.30_data=diff_processed.tif"   --focus_area 516 195 15

#v=7
# python find_defects_isolated_fixedR.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/2processed/run=20260513-032137_x=0_y=0_z=0_T=20_v=58.30_data=diff_processed.tif"  --focus_area 516 195 15


#For finding images lum vs temperature

#python find_defects_isolated_fixedR.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_v=5_run=20260514-062228/2processed/run=20260514-062228_x=0_y=0_z=0_T=17_v=56.13_data=diff_processed.tif" --focus_area 418 127 10
#python find_defects_isolated_fixedR.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_v=5_run=20260514-062228/2processed/run=20260514-062228_x=0_y=0_z=0_T=19_v=56.24_data=diff_processed.tif" --focus_area 418 127 10
#python find_defects_isolated_fixedR.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_v=5_run=20260514-062228/2processed/run=20260514-062228_x=0_y=0_z=0_T=21_v=56.35_data=diff_processed.tif" --focus_area 418 127 10
#python find_defects_isolated_fixedR.py --input "/Users/ludovicarainero/emmi_fixed_radius/DATA_irradiated/before_annealing/A1_v=5_run=20260514-062228/2processed/run=20260514-062228_x=0_y=0_z=0_T=23_v=56.46_data=diff_processed.tif" --focus_area 418 127 10


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
from photutils.segmentation import (
    SourceCatalog,
    SourceFinder,
    make_2dgaussian_kernel,
)


# ============================================================
# USER CONFIGURATION
# ============================================================

# Fixed radius used later to integrate the luminosity around
# each selected isolated source.
#
# IMPORTANT:
# This radius is NOT used for source detection and does NOT
# influence which sources are considered isolated.
integration_radius = 10.0  # pixels


# Minimum required centroid distance:
#
#   1. from every other detected source;
#   2. from every image border.
#
# A source must satisfy BOTH conditions.
isolation_radius = 60.0  # pixels


# Fixed display scale used by --focus and --focus_area.
vmin_global = 0.0
vmax_global = 5.0


# Background-estimation parameters.
background_box_size = (50, 50)
background_filter_size = (3, 3)


# Detection threshold:
#
#     threshold = threshold_sigma * local background RMS
threshold_sigma = 4.0


# Gaussian convolution kernel used for source detection.
kernel_fwhm = 3.0
kernel_size = 5


# ============================================================
# ARGUMENTS
# ============================================================

def parse_arguments():

    parser = argparse.ArgumentParser(
        description=(
            "Detect and deblend luminous sources, then select only "
            "sources isolated both from every other detected source "
            "and from the image borders."
        )
    )

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input 2D TIF filename",
    )

    parser.add_argument(
        "--output_origin",
        type=str,
        required=False,
        help="Optional PNG filename for the original image",
    )

    parser.add_argument(
        "--output",
        type=str,
        required=False,
        help=(
            "Optional PNG filename containing the processed image "
            "and the segmentation map"
        ),
    )

    parser.add_argument(
        "--display_original",
        action="store_true",
        help="Display the original image",
    )

    parser.add_argument(
        "--display_processed",
        action="store_true",
        help="Display the processed image and segmentation map",
    )

    parser.add_argument(
        "--coordinates_python",
        type=str,
        required=False,
        help=(
            "Output TXT filename for isolated-source coordinates. "
            "Format: defect_id, x, y, integration_area_pix2, "
            "integration_radius_pix"
        ),
    )

    parser.add_argument(
        "--coordinates_root",
        type=str,
        required=False,
        help=(
            "Output TXT filename for isolated-source ROOT coordinates. "
            "Format: x, y_root, integration_area_pix2, "
            "integration_radius_pix"
        ),
    )

    parser.add_argument(
        "--focus",
        type=int,
        required=False,
        help=(
            "Display a region around an isolated defect_id written "
            "in --coordinates_python"
        ),
    )

    parser.add_argument(
        "--focus_area",
        nargs=3,
        type=float,
        metavar=("X", "Y", "RADIUS"),
        help="Focus on a manually selected area: x y radius",
    )

    parser.add_argument(
        "--circles",
        action="store_true",
        help=(
            "Display all detected centroids and draw fixed "
            "integration-radius circles around isolated sources"
        ),
    )

    # --------------------------------------------------------
    # SourceFinder / deblending parameters
    # --------------------------------------------------------

    parser.add_argument(
        "--npixels",
        type=int,
        default=30,
        help="Minimum number of connected pixels for source detection",
    )

    parser.add_argument(
        "--nlevels",
        type=int,
        default=32,
        help="Number of multi-thresholding levels used for deblending",
    )

    parser.add_argument(
        "--contrast",
        type=float,
        default=0.001,
        help="Deblending contrast parameter",
    )

    return parser.parse_args()


# ============================================================
# CONFIGURATION VALIDATION
# ============================================================

def validate_configuration(args):
    """
    Validate the fixed radii and SourceFinder parameters.
    """

    if integration_radius <= 0:
        raise ValueError(
            f"Invalid integration_radius={integration_radius}. "
            "integration_radius must be greater than zero."
        )

    if isolation_radius <= 0:
        raise ValueError(
            f"Invalid isolation_radius={isolation_radius}. "
            "isolation_radius must be greater than zero."
        )

    if isolation_radius <= 2.0 * integration_radius:
        raise ValueError(
            "Invalid radius configuration: isolation_radius must be "
            "greater than 2 * integration_radius. "
            f"Current values: "
            f"isolation_radius={isolation_radius}, "
            f"integration_radius={integration_radius}, "
            f"2*integration_radius={2.0 * integration_radius}."
        )

    if args.npixels <= 0:
        raise ValueError(
            f"--npixels must be greater than zero; "
            f"got {args.npixels}."
        )

    if args.nlevels <= 0:
        raise ValueError(
            f"--nlevels must be greater than zero; "
            f"got {args.nlevels}."
        )

    if not 0.0 <= args.contrast <= 1.0:
        raise ValueError(
            f"--contrast must be between 0 and 1; "
            f"got {args.contrast}."
        )


# ============================================================
# GENERAL HELPER FUNCTIONS
# ============================================================

def get_centroid_column_names(table):
    """
    Return centroid column names for different Photutils versions.
    """

    if (
        "x_centroid" in table.colnames
        and "y_centroid" in table.colnames
    ):
        return "x_centroid", "y_centroid"

    if (
        "xcentroid" in table.colnames
        and "ycentroid" in table.colnames
    ):
        return "xcentroid", "ycentroid"

    raise RuntimeError(
        "Centroid columns were not found in the SourceCatalog table."
    )


def to_float(value):
    """
    Convert an Astropy Quantity or table value to a plain float.
    """

    if hasattr(value, "value"):
        value = value.value

    return float(value)


# ============================================================
# ISOLATED-SOURCE SELECTION
# ============================================================

def select_isolated_sources(
    table,
    x_col,
    y_col,
    image_shape,
):
    """
    Select sources isolated both from other sources and from
    the image borders.

    A source is accepted only when:

        nearest_source_distance > isolation_radius

    AND

        nearest_border_distance > isolation_radius


    Parameters
    ----------
    table : astropy.table.Table
        Complete detected/deblended source table.

    x_col, y_col : str
        Names of the centroid columns.

    image_shape : tuple
        Shape of the image as (ny, nx).


    Returns
    -------
    isolated_table : astropy.table.Table
        Table containing only isolated sources.

    nearest_source_distance : numpy.ndarray
        Distance from the nearest detected source for every source.

    border_distance : numpy.ndarray
        Minimum distance from the four image borders for every source.

    isolated_mask : numpy.ndarray
        Boolean mask relative to the complete detected-source table.
    """

    ny, nx = image_shape

    # --------------------------------------------------------
    # Extract centroid coordinates
    # --------------------------------------------------------

    x = np.array(
        [to_float(value) for value in table[x_col]],
        dtype=float,
    )

    y = np.array(
        [to_float(value) for value in table[y_col]],
        dtype=float,
    )

    number_of_sources = len(table)

    # ========================================================
    # 1. DISTANCE FROM OTHER SOURCES
    # ========================================================

    if number_of_sources == 1:

        # If there is only one source, there are no neighbours.
        nearest_source_distance = np.array(
            [np.inf],
            dtype=float,
        )

    else:

        coordinates = np.column_stack((x, y))

        # differences[i,j] contains:
        #
        #     position_i - position_j
        #
        differences = (
            coordinates[:, np.newaxis, :]
            - coordinates[np.newaxis, :, :]
        )

        # Euclidean distance between every pair of centroids.
        distance_matrix = np.sqrt(
            np.sum(differences**2, axis=2)
        )

        # Ignore each source's distance from itself.
        np.fill_diagonal(
            distance_matrix,
            np.inf,
        )

        # Minimum centroid-to-centroid distance for every source.
        nearest_source_distance = np.min(
            distance_matrix,
            axis=1,
        )

    # Source must be farther than isolation_radius
    # from every other source.
    isolated_from_sources = (
        nearest_source_distance > isolation_radius
    )

    # ========================================================
    # 2. DISTANCE FROM IMAGE BORDERS
    # ========================================================

    # Pixel coordinates range from:
    #
    #     x = 0 ... nx - 1
    #     y = 0 ... ny - 1

    distance_left = x

    distance_right = (
        (nx - 1) - x
    )

    distance_top = y

    distance_bottom = (
        (ny - 1) - y
    )

    # Minimum distance from any of the four borders.
    border_distance = np.minimum.reduce(
        [
            distance_left,
            distance_right,
            distance_top,
            distance_bottom,
        ]
    )

    # Source must also be farther than isolation_radius
    # from every border.
    isolated_from_border = (
        border_distance > isolation_radius
    )

    # ========================================================
    # FINAL ISOLATION CONDITION
    # ========================================================

    isolated_mask = (
        isolated_from_sources
        & isolated_from_border
    )

    isolated_table = table[isolated_mask].copy()

    # --------------------------------------------------------
    # Add diagnostic information to isolated-source table
    # --------------------------------------------------------

    isolated_table[
        "nearest_source_distance_pix"
    ] = nearest_source_distance[isolated_mask]

    isolated_table[
        "nearest_border_distance_pix"
    ] = border_distance[isolated_mask]

    # Consecutive defect IDs referring ONLY to isolated sources.
    isolated_table["defect_id"] = np.arange(
        1,
        len(isolated_table) + 1,
    )

    return (
        isolated_table,
        nearest_source_distance,
        border_distance,
        isolated_mask,
    )


# ============================================================
# FOCUS DISPLAY
# ============================================================

def display_focus(
    processed,
    xc,
    yc,
    radius,
    title,
):
    """
    Display a cropped image around a selected point.
    """

    half_size = max(
        20,
        int(np.ceil(2.0 * radius)),
    )

    ny, nx = processed.shape

    x_min = max(
        0,
        int(round(xc)) - half_size,
    )

    x_max = min(
        nx,
        int(round(xc)) + half_size,
    )

    y_min = max(
        0,
        int(round(yc)) - half_size,
    )

    y_max = min(
        ny,
        int(round(yc)) + half_size,
    )

    focus_image = processed[
        y_min:y_max,
        x_min:x_max,
    ]

    if focus_image.size == 0:
        raise ValueError(
            "The requested focus area around "
            f"x={xc}, y={yc} is outside the image."
        )

    # Coordinates relative to the cropped image.
    xc_local = xc - x_min
    yc_local = yc - y_min

    norm = ImageNormalize(
        vmin=vmin_global,
        vmax=vmax_global,
        stretch=SqrtStretch(),
    )

    fig_focus, ax_focus = plt.subplots(
        figsize=(6, 6)
    )

    ax_focus.imshow(
        focus_image,
        origin="upper",
        norm=norm,
    )

    circle = Circle(
        (xc_local, yc_local),
        radius=radius,
        edgecolor="red",
        facecolor="none",
        linewidth=1.5,
    )

    ax_focus.add_patch(circle)

    ax_focus.set_title(title)

    ax_focus.set_xlim(
        0,
        focus_image.shape[1],
    )

    ax_focus.set_ylim(
        focus_image.shape[0],
        0,
    )

    plt.tight_layout()
    plt.show()
    plt.close(fig_focus)


# ============================================================
# ORIGINAL IMAGE DISPLAY / SAVE
# ============================================================

def save_original_image(
    image,
    args,
):
    """
    Display and/or save the original image.
    """

    if not (
        args.display_original
        or args.output_origin
    ):
        return

    fig, ax = plt.subplots(
        figsize=(10, 5)
    )

    ax.imshow(
        image,
        origin="upper",
    )

    ax.axis("off")

    if args.output_origin:

        fig.savefig(
            args.output_origin,
            format="png",
            dpi=300,
            bbox_inches="tight",
            pad_inches=0,
        )

    if args.display_original:
        plt.show()

    plt.close(fig)


# ============================================================
# PROCESSED IMAGE DISPLAY / SAVE
# ============================================================

def save_processed_images(
    processed,
    segment_map,
    args,
):
    """
    Display and/or save the background-subtracted image and
    the complete deblended segmentation map.
    """

    if not (
        args.display_processed
        or args.output
    ):
        return

    norm = ImageNormalize(
        stretch=SqrtStretch()
    )

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=(10, 12.5),
    )

    # --------------------------------------------------------
    # Background-subtracted image
    # --------------------------------------------------------

    ax1.imshow(
        processed,
        origin="upper",
        norm=norm,
    )

    ax1.set_title(
        "Background-subtracted Data"
    )

    ax1.axis("off")

    # --------------------------------------------------------
    # Segmentation map
    # --------------------------------------------------------

    ax2.imshow(
        segment_map.data,
        origin="upper",
        cmap=segment_map.cmap,
        interpolation="nearest",
    )

    ax2.set_title(
        "Deblended Segmentation Image — All Detected Sources"
    )

    ax2.axis("off")

    plt.tight_layout()

    if args.output:

        fig.savefig(
            args.output,
            dpi=300,
            bbox_inches="tight",
        )

    if args.display_processed:
        plt.show()

    plt.close(fig)


# ============================================================
# COORDINATE OUTPUT
# ============================================================

def save_python_coordinates(
    isolated_table,
    x_col,
    y_col,
    filename,
):
    """
    Save isolated-source coordinates for Python.

    Format:

        defect_id,
        x,
        y,
        integration_area_pix2,
        integration_radius_pix
    """

    fixed_area = (
        np.pi * integration_radius**2
    )

    print(
        " --- saving isolated Python coordinates to:",
        filename,
    )

    with open(
        filename,
        "w",
        encoding="utf-8",
    ) as output_file:

        for row in isolated_table:

            defect_id = int(
                row["defect_id"]
            )

            x = to_float(
                row[x_col]
            )

            y = to_float(
                row[y_col]
            )

            output_file.write(
                f"{defect_id}, "
                f"{x:.2f}, "
                f"{y:.2f}, "
                f"{fixed_area:.2f}, "
                f"{integration_radius:.2f}\n"
            )


def save_root_coordinates(
    isolated_table,
    x_col,
    y_col,
    image_height,
    filename,
):
    """
    Save isolated-source coordinates for ROOT.

    Format:

        x,
        y_root,
        integration_area_pix2,
        integration_radius_pix

    where:

        y_root = image_height - y_python
    """

    fixed_area = (
        np.pi * integration_radius**2
    )

    print(
        " --- saving isolated ROOT coordinates to:",
        filename,
    )

    with open(
        filename,
        "w",
        encoding="utf-8",
    ) as output_file:

        for row in isolated_table:

            x = to_float(
                row[x_col]
            )

            y_python = to_float(
                row[y_col]
            )

            y_root = (
                float(image_height)
                - y_python
            )

            output_file.write(
                f"{x:.2f}, "
                f"{y_root:.2f}, "
                f"{fixed_area:.2f}, "
                f"{integration_radius:.2f}\n"
            )


# ============================================================
# CIRCLE DISPLAY
# ============================================================

def display_source_circles(
    processed,
    all_sources,
    isolated_sources,
    x_col,
    y_col,
):
    """
    Display:

    - all detected centroids as yellow crosses;
    - fixed integration-radius circles around isolated sources.
    """

    norm = ImageNormalize(
        stretch=SqrtStretch()
    )

    fig, ax = plt.subplots(
        figsize=(8, 8)
    )

    ax.imshow(
        processed,
        origin="upper",
        norm=norm,
    )

    ax.set_title(
        "Detected Sources and Isolated Integration Regions\n"
        f"integration radius = {integration_radius:.1f} px, "
        f"isolation radius = {isolation_radius:.1f} px"
    )

    # --------------------------------------------------------
    # All detected centroids
    # --------------------------------------------------------

    for row in all_sources:

        x = to_float(
            row[x_col]
        )

        y = to_float(
            row[y_col]
        )

        ax.plot(
            x,
            y,
            marker="x",
            markersize=4,
            color="yellow",
        )

    # --------------------------------------------------------
    # Isolated sources
    # --------------------------------------------------------

    for row in isolated_sources:

        defect_id = int(
            row["defect_id"]
        )

        x = to_float(
            row[x_col]
        )

        y = to_float(
            row[y_col]
        )

        circle = Circle(
            (x, y),
            radius=integration_radius,
            edgecolor="red",
            facecolor="none",
            linewidth=1.5,
        )

        ax.add_patch(circle)

        ax.text(
            x,
            y,
            str(defect_id),
            color="red",
            fontsize=8,
        )

    ax.set_xlim(
        0,
        processed.shape[1],
    )

    ax.set_ylim(
        processed.shape[0],
        0,
    )

    plt.tight_layout()
    plt.show()
    plt.close(fig)


# ============================================================
# MAIN PROGRAM
# ============================================================

def main():

    args = parse_arguments()

    # ========================================================
    # Validate configuration
    # ========================================================

    try:

        validate_configuration(args)

    except ValueError as error:

        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )

        return 1


    fixed_area = (
        np.pi * integration_radius**2
    )

    print(
        f" --- integration_radius = "
        f"{integration_radius:.2f} pixels"
    )

    print(
        f" --- integration area   = "
        f"{fixed_area:.2f} pixels^2"
    )

    print(
        f" --- isolation_radius   = "
        f"{isolation_radius:.2f} pixels"
    )


    # ========================================================
    # Read TIF image
    # ========================================================

    print(
        " --- opening input image:",
        args.input,
    )

    try:

        image = tifffile.imread(
            args.input
        ).astype(float)

    except (
        OSError,
        ValueError,
        tifffile.TiffFileError,
    ) as error:

        print(
            f"ERROR: could not read input image: {error}",
            file=sys.stderr,
        )

        return 1


    if image.ndim != 2:

        print(
            "ERROR: expected a 2D image, "
            f"but input shape is {image.shape}.",
            file=sys.stderr,
        )

        return 1


    # ========================================================
    # Original image
    # ========================================================

    save_original_image(
        image,
        args,
    )


    # ========================================================
    # Background subtraction
    # ========================================================

    background_estimator = (
        MedianBackground()
    )

    try:

        background = Background2D(
            image,
            box_size=background_box_size,
            filter_size=background_filter_size,
            bkg_estimator=background_estimator,
        )

    except ValueError as error:

        print(
            f"ERROR: background estimation failed: {error}",
            file=sys.stderr,
        )

        return 1


    processed = (
        image
        - background.background
    )


    # ========================================================
    # Detection threshold
    # ========================================================

    threshold = (
        threshold_sigma
        * background.background_rms
    )


    # ========================================================
    # Gaussian convolution
    # ========================================================

    kernel = make_2dgaussian_kernel(
        kernel_fwhm,
        size=kernel_size,
    )

    convolved_data = convolve(
        processed,
        kernel,
    )


    # ========================================================
    # Detection and deblending
    # ========================================================

    finder = SourceFinder(
        npixels=args.npixels,
        deblend=True,
        nlevels=args.nlevels,
        contrast=args.contrast,
        progress_bar=False,
    )

    segment_map = finder(
        convolved_data,
        threshold,
    )


    # ========================================================
    # No sources detected
    # ========================================================

    if segment_map is None:

        print(
            " --- No sources detected."
        )

        # Create empty coordinate files if requested,
        # avoiding obsolete files from previous analyses.

        if args.coordinates_python:

            open(
                args.coordinates_python,
                "w",
                encoding="utf-8",
            ).close()

            print(
                " --- empty Python coordinate file created:",
                args.coordinates_python,
            )


        if args.coordinates_root:

            open(
                args.coordinates_root,
                "w",
                encoding="utf-8",
            ).close()

            print(
                " --- empty ROOT coordinate file created:",
                args.coordinates_root,
            )

        return 0


    # ========================================================
    # Source catalogue
    # ========================================================

    catalog = SourceCatalog(
        processed,
        segment_map,
        convolved_data=convolved_data,
    )

    table = catalog.to_table()

    x_col, y_col = (
        get_centroid_column_names(table)
    )

    table[x_col].info.format = ".2f"
    table[y_col].info.format = ".2f"

    # Consecutive ID for the COMPLETE detected-source list.
    table["all_source_id"] = np.arange(
        1,
        len(table) + 1,
    )


    # ========================================================
    # Select isolated sources
    # ========================================================

    (
        isolated_table,
        nearest_source_distances,
        border_distances,
        isolated_mask,
    ) = select_isolated_sources(
        table,
        x_col,
        y_col,
        processed.shape,
    )


    # ========================================================
    # Statistics
    # ========================================================

    number_all = len(table)

    number_isolated = len(
        isolated_table
    )

    number_rejected_neighbours = int(
        np.count_nonzero(
            nearest_source_distances
            <= isolation_radius
        )
    )

    number_rejected_border = int(
        np.count_nonzero(
            border_distances
            <= isolation_radius
        )
    )


    print(
        f" --- detected and deblended sources: "
        f"{number_all}"
    )

    print(
        f" --- isolated sources:              "
        f"{number_isolated}"
    )

    print(
        " --- sources failing neighbour condition: "
        f"{number_rejected_neighbours}"
    )

    print(
        " --- sources failing border condition:    "
        f"{number_rejected_border}"
    )


    # ========================================================
    # Print isolated-source information
    # ========================================================

    if number_isolated > 0:

        print(
            " --- isolated-source centroids:"
        )

        for row in isolated_table:

            print(
                f"     defect_id="
                f"{int(row['defect_id'])}, "
                f"x={to_float(row[x_col]):.2f}, "
                f"y={to_float(row[y_col]):.2f}, "
                f"nearest source="
                f"{to_float(row['nearest_source_distance_pix']):.2f} px, "
                f"nearest border="
                f"{to_float(row['nearest_border_distance_pix']):.2f} px"
            )


    # ========================================================
    # Processed image / segmentation map
    # ========================================================

    save_processed_images(
        processed,
        segment_map,
        args,
    )


    # ========================================================
    # Focus on isolated source
    # ========================================================

    if args.focus is not None:

        matches = np.where(
            np.array(
                isolated_table["defect_id"],
                dtype=int,
            )
            == args.focus
        )[0]


        if len(matches) == 0:

            if number_isolated == 0:

                print(
                    "ERROR: no isolated sources are available "
                    "for --focus."
                )

            else:

                print(
                    f"ERROR: selected isolated defect_id "
                    f"{args.focus} is out of range. "
                    f"Available values: 1 ... "
                    f"{number_isolated}"
                )

            return 1


        focus_index = int(
            matches[0]
        )

        xc = to_float(
            isolated_table[x_col][focus_index]
        )

        yc = to_float(
            isolated_table[y_col][focus_index]
        )

        nearest_source = to_float(
            isolated_table[
                "nearest_source_distance_pix"
            ][focus_index]
        )

        nearest_border = to_float(
            isolated_table[
                "nearest_border_distance_pix"
            ][focus_index]
        )


        print(
            f" --- focus on isolated defect "
            f"{args.focus}: "
            f"x={xc:.2f}, "
            f"y={yc:.2f}, "
            f"integration radius="
            f"{integration_radius:.2f} px, "
            f"nearest source="
            f"{nearest_source:.2f} px, "
            f"nearest border="
            f"{nearest_border:.2f} px"
        )


        display_focus(
            processed,
            xc,
            yc,
            integration_radius,
            f"Isolated defect 15 with radius {integration_radius:.1f} px",
        )


    # ========================================================
    # Focus on manually selected area
    # ========================================================

    if args.focus_area is not None:

        xc, yc, manual_radius = (
            args.focus_area
        )

        if manual_radius <= 0:

            print(
                "ERROR: the --focus_area radius "
                "must be greater than zero."
            )

            return 1


        print(
            f" --- focus on manual area: "
            f"x={xc:.2f}, "
            f"y={yc:.2f}, "
            f"radius={manual_radius:.2f} pixels"
        )


        display_focus(
            processed,
            xc,
            yc,
            manual_radius,
            (
                f" Hotspot 25: "
                f"x={xc:.1f}, "
                f"y={yc:.1f}, "
                f"r={manual_radius:.1f}px, "
                f" T=23.0°C "
            ),
        )


    # ========================================================
    # Save coordinates
    # ========================================================

    if args.coordinates_python:

        save_python_coordinates(
            isolated_table,
            x_col,
            y_col,
            args.coordinates_python,
        )


    if args.coordinates_root:

        save_root_coordinates(
            isolated_table,
            x_col,
            y_col,
            processed.shape[0],
            args.coordinates_root,
        )


    # ========================================================
    # Display circles
    # ========================================================

    if args.circles:

        display_source_circles(
            processed,
            table,
            isolated_table,
            x_col,
            y_col,
        )


    return 0


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    sys.exit(main())