#!/usr/bin/env python3
"""
Detect luminous sources in a 2D TIF image using Photutils.

Processing steps
----------------
1. Estimate and subtract the median background with Background2D.
2. Convolve the background-subtracted image with a 2D Gaussian kernel.
3. Detect and deblend sources with SourceFinder.
4. Calculate each source area from the segmentation footprint.
5. Calculate a variable equivalent radius:

       equivalent_radius = sqrt(segment_area / pi)

6. Select only isolated sources. A source is isolated when its centroid is:
   - farther than isolation_radius from every other detected centroid;
   - farther than isolation_radius from each of the four image borders.
7. With --coordinates_root, save only the isolated sources as:

       x, y_root, segment_area_pix2, equivalent_radius_pix

   where y_root = image_height - y_python.
"""

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import tifffile

from astropy.convolution import convolve
from astropy.visualization import ImageNormalize, SqrtStretch
from matplotlib.patches import Circle
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import SourceCatalog, SourceFinder, make_2dgaussian_kernel


# ============================================================
# USER CONFIGURATION
# ============================================================
# Minimum required centroid distance both from every other detected
# source and from every image border.
isolation_radius = 60.0  # pixels

# Background-estimation parameters.
background_box_size = (50, 50)
background_filter_size = (3, 3)

# Detection threshold = threshold_sigma * local background RMS.
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
            "Detect and deblend luminous sources, calculate a variable source "
            "radius from the segmented area, and save only isolated sources."
        )
    )

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input 2D TIF filename",
    )
    parser.add_argument(
        "--coordinates_root",
        type=str,
        required=False,
        help=(
            "Output TXT filename. Only isolated sources are written as: "
            "x, y_root, segment_area_pix2, equivalent_radius_pix"
        ),
    )

    parser.add_argument(
        "--output_origin",
        type=str,
        required=False,
        help="Optional output PNG filename for the original image",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=False,
        help="Optional output PNG filename for the processed and segmented images",
    )
    parser.add_argument(
        "--display_original",
        action="store_true",
        help="Display the original image",
    )
    parser.add_argument(
        "--display_processed",
        action="store_true",
        help="Display the background-subtracted image and segmentation map",
    )
    parser.add_argument(
        "--circles",
        action="store_true",
        help=(
            "Display all detected centroids and variable-radius circles around "
            "the isolated sources"
        ),
    )

    # SourceFinder / deblending parameters.
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
# HELPER FUNCTIONS
# ============================================================
def validate_configuration(args):
    """Check configuration and command-line values before processing."""
    if isolation_radius <= 0:
        raise ValueError(
            f"isolation_radius must be greater than zero; got {isolation_radius}."
        )

    if args.npixels <= 0:
        raise ValueError(f"--npixels must be greater than zero; got {args.npixels}.")

    if args.nlevels <= 0:
        raise ValueError(f"--nlevels must be greater than zero; got {args.nlevels}.")

    if not 0.0 <= args.contrast <= 1.0:
        raise ValueError(
            f"--contrast must be between 0 and 1; got {args.contrast}."
        )


def to_float(value):
    """Convert an Astropy Quantity or table value to a plain float."""
    if hasattr(value, "value"):
        value = value.value
    return float(value)


def get_centroid_column_names(table):
    """Support both current and older Photutils centroid column names."""
    if "x_centroid" in table.colnames and "y_centroid" in table.colnames:
        return "x_centroid", "y_centroid"

    if "xcentroid" in table.colnames and "ycentroid" in table.colnames:
        return "xcentroid", "ycentroid"

    raise RuntimeError("Centroid columns were not found in the SourceCatalog table.")


def calculate_source_dimensions(catalog, table):
    """
    Add segmented source area and equivalent circular radius to the table.

    The equivalent radius is calculated independently for every source as:

        radius = sqrt(segment_area / pi)
    """
    segment_area = np.array(
        [to_float(area) for area in catalog.segment_area],
        dtype=float,
    )

    equivalent_radius = np.sqrt(segment_area / np.pi)

    table["segment_area_pix2"] = segment_area
    table["equivalent_radius_pix"] = equivalent_radius

    table["segment_area_pix2"].info.format = ".2f"
    table["equivalent_radius_pix"].info.format = ".2f"


def select_isolated_sources(table, x_col, y_col, image_shape):
    """
    Select sources isolated from both other sources and the image borders.

    The distance from other sources is the Euclidean centroid-to-centroid
    distance. The border distance is the minimum centroid distance from the
    left, right, top, and bottom borders.

    Returns
    -------
    isolated_table : astropy.table.Table
        Table containing only isolated sources.
    isolated_mask : numpy.ndarray
        Boolean mask relative to the complete source table.
    nearest_source_distance : numpy.ndarray
        Nearest-neighbour distance for each detected source.
    border_distance : numpy.ndarray
        Minimum border distance for each detected source.
    """
    ny, nx = image_shape

    x = np.array([to_float(value) for value in table[x_col]], dtype=float)
    y = np.array([to_float(value) for value in table[y_col]], dtype=float)

    number_of_sources = len(table)

    # --------------------------------------------------------
    # Distance from every other detected source
    # --------------------------------------------------------
    if number_of_sources == 1:
        nearest_source_distance = np.array([np.inf], dtype=float)
    else:
        coordinates = np.column_stack((x, y))
        differences = (
            coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
        )
        distance_matrix = np.sqrt(np.sum(differences**2, axis=2))

        # Ignore the zero distance of each source from itself.
        np.fill_diagonal(distance_matrix, np.inf)
        nearest_source_distance = np.min(distance_matrix, axis=1)

    isolated_from_sources = nearest_source_distance > isolation_radius

    # --------------------------------------------------------
    # Distance from the four image borders
    # --------------------------------------------------------
    # Pixel coordinates extend from 0 to nx - 1 and from 0 to ny - 1.
    distance_left = x
    distance_right = (nx - 1) - x
    distance_top = y
    distance_bottom = (ny - 1) - y

    border_distance = np.minimum.reduce(
        [distance_left, distance_right, distance_top, distance_bottom]
    )
    isolated_from_border = border_distance > isolation_radius

    # A source must pass both isolation conditions.
    isolated_mask = isolated_from_sources & isolated_from_border

    isolated_table = table[isolated_mask].copy()
    isolated_table["nearest_source_distance_pix"] = nearest_source_distance[
        isolated_mask
    ]
    isolated_table["nearest_border_distance_pix"] = border_distance[isolated_mask]

    # Consecutive ID referring only to the final isolated-source list.
    isolated_table["isolated_source_id"] = np.arange(
        1,
        len(isolated_table) + 1,
    )

    return (
        isolated_table,
        isolated_mask,
        nearest_source_distance,
        border_distance,
    )


def save_original_image(image, args):
    """Display and/or save the original image."""
    if not (args.display_original or args.output_origin):
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.imshow(image, origin="upper")
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


def save_processed_images(processed, segment_map, args):
    """Display and/or save the processed image and deblended segmentation map."""
    if not (args.display_processed or args.output):
        return

    norm = ImageNormalize(stretch=SqrtStretch())

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))

    ax1.imshow(processed, origin="upper", norm=norm)
    ax1.set_title("Background-subtracted Data")
    ax1.axis("off")

    ax2.imshow(
        segment_map.data,
        origin="upper",
        cmap=segment_map.cmap,
        interpolation="nearest",
    )
    ax2.set_title("Deblended Segmentation Image — All Detected Sources")
    ax2.axis("off")

    plt.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=300, bbox_inches="tight")

    if args.display_processed:
        plt.show()

    plt.close(fig)


def display_source_circles(processed, all_sources, isolated_sources, x_col, y_col):
    """Display all centroids and variable-radius circles around isolated sources."""
    norm = ImageNormalize(stretch=SqrtStretch())

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(processed, origin="upper", norm=norm)
    ax.set_title(
        "Detected Sources and Isolated Variable-Radius Regions\n"
        f"isolation radius = {isolation_radius:.1f} px"
    )

    # Mark every detected and deblended source centroid.
    for row in all_sources:
        x = to_float(row[x_col])
        y = to_float(row[y_col])
        ax.plot(x, y, marker="x", markersize=4, color="yellow")

    # Draw each isolated source using its own equivalent radius.
    for row in isolated_sources:
        source_id = int(row["isolated_source_id"])
        x = to_float(row[x_col])
        y = to_float(row[y_col])
        radius = to_float(row["equivalent_radius_pix"])

        circle = Circle(
            (x, y),
            radius=radius,
            edgecolor="red",
            facecolor="none",
            linewidth=1.5,
        )
        ax.add_patch(circle)
        ax.text(x, y, str(source_id), color="red", fontsize=8)

    ax.set_xlim(0, processed.shape[1])
    ax.set_ylim(processed.shape[0], 0)

    plt.tight_layout()
    plt.show()
    plt.close(fig)


def save_root_coordinates(isolated_sources, x_col, y_col, image_height, filename):
    """
    Save only isolated sources in the coordinate format expected by ROOT.

    Columns:
        x, y_root, segment_area_pix2, equivalent_radius_pix
    """
    print(" --- saving isolated ROOT coordinates to:", filename)

    with open(filename, "w", encoding="utf-8") as output_file:
        for row in isolated_sources:
            x = to_float(row[x_col])
            y_python = to_float(row[y_col])
            y_root = float(image_height) - y_python
            area = to_float(row["segment_area_pix2"])
            radius = to_float(row["equivalent_radius_pix"])

            output_file.write(
                f"{x:.2f}, {y_root:.2f}, {area:.2f}, {radius:.2f}\n"
            )


# ============================================================
# MAIN PROGRAM
# ============================================================
def main():
    args = parse_arguments()

    try:
        validate_configuration(args)
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    # ------------------------------------------------------------
    # Read the TIF image
    # ------------------------------------------------------------
    print(" --- opening input image:", args.input)

    try:
        image = tifffile.imread(args.input).astype(float)
    except (OSError, ValueError, tifffile.TiffFileError) as error:
        print(f"ERROR: could not read input image: {error}", file=sys.stderr)
        return 1

    if image.ndim != 2:
        print(
            f"ERROR: expected a 2D image, but the input shape is {image.shape}.",
            file=sys.stderr,
        )
        return 1

    save_original_image(image, args)

    # ------------------------------------------------------------
    # Median-background estimation and subtraction
    # ------------------------------------------------------------
    background_estimator = MedianBackground()

    try:
        background = Background2D(
            image,
            box_size=background_box_size,
            filter_size=background_filter_size,
            bkg_estimator=background_estimator,
        )
    except ValueError as error:
        print(f"ERROR: background estimation failed: {error}", file=sys.stderr)
        return 1

    processed = image - background.background

    # Local detection threshold, as in the first script.
    threshold = threshold_sigma * background.background_rms

    # ------------------------------------------------------------
    # Gaussian convolution — always performed and used for detection
    # ------------------------------------------------------------
    kernel = make_2dgaussian_kernel(kernel_fwhm, size=kernel_size)
    convolved_data = convolve(processed, kernel)

    # ------------------------------------------------------------
    # Segmentation and deblending of all sources
    # ------------------------------------------------------------
    finder = SourceFinder(
        npixels=args.npixels,
        deblend=True,
        nlevels=args.nlevels,
        contrast=args.contrast,
        progress_bar=False,
    )

    segment_map = finder(convolved_data, threshold)

    if segment_map is None:
        print(" --- No sources detected.")

        # Create an empty requested coordinate file rather than leaving an
        # obsolete result from an earlier run in place.
        if args.coordinates_root:
            open(args.coordinates_root, "w", encoding="utf-8").close()
            print(" --- empty coordinate file created:", args.coordinates_root)

        return 0

    # ------------------------------------------------------------
    # Source catalogue and variable source dimensions
    # ------------------------------------------------------------
    catalog = SourceCatalog(
        processed,
        segment_map,
        convolved_data=convolved_data,
    )

    source_table = catalog.to_table()
    x_col, y_col = get_centroid_column_names(source_table)

    source_table[x_col].info.format = ".2f"
    source_table[y_col].info.format = ".2f"
    source_table["all_source_id"] = np.arange(1, len(source_table) + 1)

    calculate_source_dimensions(catalog, source_table)

    # ------------------------------------------------------------
    # Select sources isolated from neighbours and image borders
    # ------------------------------------------------------------
    (
        isolated_table,
        isolated_mask,
        nearest_source_distances,
        border_distances,
    ) = select_isolated_sources(
        source_table,
        x_col,
        y_col,
        processed.shape,
    )

    number_all = len(source_table)
    number_isolated = len(isolated_table)
    number_rejected_neighbours = int(
        np.count_nonzero(nearest_source_distances <= isolation_radius)
    )
    number_rejected_border = int(
        np.count_nonzero(border_distances <= isolation_radius)
    )

    print(f" --- detected and deblended sources: {number_all}")
    print(f" --- isolated sources:              {number_isolated}")
    print(
        " --- sources failing neighbour condition: "
        f"{number_rejected_neighbours}"
    )
    print(
        " --- sources failing border condition:    "
        f"{number_rejected_border}"
    )

    if number_isolated > 0:
        print(" --- isolated sources:")
        for row in isolated_table:
            print(
                f"     id={int(row['isolated_source_id'])}, "
                f"x={to_float(row[x_col]):.2f}, "
                f"y={to_float(row[y_col]):.2f}, "
                f"area={to_float(row['segment_area_pix2']):.2f} pix^2, "
                f"radius={to_float(row['equivalent_radius_pix']):.2f} pix, "
                "nearest source="
                f"{to_float(row['nearest_source_distance_pix']):.2f} pix, "
                "nearest border="
                f"{to_float(row['nearest_border_distance_pix']):.2f} pix"
            )

    # ------------------------------------------------------------
    # Optional outputs
    # ------------------------------------------------------------
    save_processed_images(processed, segment_map, args)

    if args.coordinates_root:
        save_root_coordinates(
            isolated_table,
            x_col,
            y_col,
            processed.shape[0],
            args.coordinates_root,
        )

    if args.circles:
        display_source_circles(
            processed,
            source_table,
            isolated_table,
            x_col,
            y_col,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
