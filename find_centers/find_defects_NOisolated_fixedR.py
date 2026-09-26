#!/usr/bin/env python3

"""
Detect luminous defects in an EMMI TIF image.

Main characteristics of this version:
  1. ALL detected/deblended sources are kept (no isolation selection).
  2. A FIXED integration radius is assigned to every source.
  3. The fixed radius is R = 20 pixels.
  4. --circle / --circles displays all detected defects with R = 20 px.
  5. --save saves figures produced by display commands in both PNG and PDF.
  6. --device adds an optional device label to figure titles.

The coordinate outputs are:

--coordinates_python
    defect_id, x, y, integration_area_pix2, integration_radius_pix

--coordinates_root
    x, y_root, integration_area_pix2, integration_radius_pix
"""

import argparse
import sys
from pathlib import Path

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

# Fixed radius assigned to EVERY detected source.
integration_radius = 20.0  # pixels

# Fixed display scale used by --focus and --focus_area.
vmin_global = 0.0
vmax_global = 5.0

# Background-estimation parameters.
background_box_size = (50, 50)
background_filter_size = (3, 3)

# Detection threshold:
# threshold = threshold_sigma * local background RMS
threshold_sigma = 3.0

# Gaussian convolution kernel used for source detection.
kernel_fwhm = 3.0
kernel_size = 5


# ============================================================
# ARGUMENTS
# ============================================================

def parse_arguments():

    parser = argparse.ArgumentParser(
        description=(
            "Detect and deblend ALL luminous sources in an EMMI image "
            "and assign the same fixed integration radius R=20 px "
            "to every detected defect."
        )
    )

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input 2D TIF filename",
    )

    parser.add_argument(
        "--device",
        type=str,
        required=False,
        default=None,
        help=(
            "Optional device label added at the beginning of figure titles. "
            "Use quotes when the label contains spaces, e.g. "
            '--device "A1 new device"'
        ),
    )

    parser.add_argument(
        "--output_origin",
        type=str,
        required=False,
        help="Optional output filename for the original image",
    )

    parser.add_argument(
        "--output",
        type=str,
        required=False,
        help=(
            "Optional output filename containing the processed image "
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

    # Preserved from find_defects_NOisolated_changeR.py.
    # If omitted, the processed (non-convolved) image is used for detection.
    parser.add_argument(
        "--convolution",
        action="store_true",
        help="Apply Gaussian convolution to the image used for source detection",
    )

    parser.add_argument(
        "--coordinates_python",
        type=str,
        required=False,
        help=(
            "Output TXT filename for ALL detected sources. "
            "Format: defect_id, x, y, integration_area_pix2, "
            "integration_radius_pix"
        ),
    )

    parser.add_argument(
        "--coordinates_root",
        type=str,
        required=False,
        help=(
            "Output TXT filename for ALL detected sources in ROOT coordinates. "
            "Format: x, y_root, integration_area_pix2, "
            "integration_radius_pix"
        ),
    )

    parser.add_argument(
        "--focus",
        type=int,
        required=False,
        help=(
            "Display a region around the selected defect_id written "
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

    # Singular and plural forms are accepted for compatibility.
    parser.add_argument(
        "--circle",
        "--circles",
        dest="circles",
        action="store_true",
        help=(
            "Display the background-subtracted image with a fixed "
            "R=20 px circle around every detected source"
        ),
    )

    parser.add_argument(
        "--save",
        action="store_true",
        help=(
            "Save each figure produced by a display command in both "
            "PNG and PDF format"
        ),
    )

    parser.add_argument(
        "--save_dir",
        type=str,
        required=False,
        default=None,
        help=(
            "Directory used by --save. If omitted, figures are saved "
            "in the same directory as the input TIF"
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

    if integration_radius <= 0:
        raise ValueError(
            f"Invalid integration_radius={integration_radius}. "
            "integration_radius must be greater than zero."
        )

    if args.npixels <= 0:
        raise ValueError(
            f"--npixels must be greater than zero; got {args.npixels}."
        )

    if args.nlevels <= 0:
        raise ValueError(
            f"--nlevels must be greater than zero; got {args.nlevels}."
        )

    if not 0.0 <= args.contrast <= 1.0:
        raise ValueError(
            f"--contrast must be between 0 and 1; got {args.contrast}."
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


def get_save_directory(args):
    """
    Return the directory used by --save.
    """

    if args.save_dir is not None:
        save_dir = Path(args.save_dir).expanduser()
    else:
        save_dir = Path(args.input).expanduser().resolve().parent

    save_dir.mkdir(parents=True, exist_ok=True)

    return save_dir


def get_input_stem(args):
    """
    Return input filename without extension.
    """

    return Path(args.input).stem


def figure_title(args, title):
    """
    Add the optional device label to a figure title.

    Example:
        --device "A1 new device"

    gives:
        A1 new device - Background subtracted data with detected sources
    """

    if args.device:
        device = args.device.strip()

        if device:
            return f"{device} - {title}"

    return title


def save_figure_png_pdf(fig, args, suffix):
    """
    Save a Matplotlib figure in both PNG and vector PDF format.

    Files are written only when --save is active.
    """

    if not args.save:
        return

    save_dir = get_save_directory(args)
    base = save_dir / f"{get_input_stem(args)}_{suffix}"

    png_filename = base.with_suffix(".png")
    pdf_filename = base.with_suffix(".pdf")

    fig.savefig(
        png_filename,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf_filename,
        bbox_inches="tight",
    )

    print(" --- saved:", png_filename)
    print(" --- saved:", pdf_filename)


def save_requested_output(fig, filename):
    """
    Preserve the legacy --output / --output_origin behavior.

    The format is inferred from the filename extension. If no extension
    is supplied, PNG is used.
    """

    if not filename:
        return

    output_path = Path(filename)

    if output_path.suffix == "":
        output_path = output_path.with_suffix(".png")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    kwargs = {
        "bbox_inches": "tight",
    }

    if output_path.suffix.lower() == ".png":
        kwargs["dpi"] = 300

    fig.savefig(
        output_path,
        **kwargs,
    )

    print(" --- saved:", output_path)


# ============================================================
# FOCUS DISPLAY / SAVE
# ============================================================

def display_focus(
    processed,
    xc,
    yc,
    radius,
    title,
    args,
    save_suffix,
):
    """
    Display a cropped image around a selected point.

    The supplied radius controls only the displayed circle/crop.
    For --focus on a detected defect, radius = integration_radius = 20 px.
    For --focus_area, the manually supplied radius is preserved.
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

    ax_focus.set_title(figure_title(args, title))

    ax_focus.set_xlim(
        0,
        focus_image.shape[1],
    )

    ax_focus.set_ylim(
        focus_image.shape[0],
        0,
    )

    plt.tight_layout()

    save_figure_png_pdf(
        fig_focus,
        args,
        save_suffix,
    )

    plt.show()
    plt.close(fig_focus)


# ============================================================
# ORIGINAL IMAGE DISPLAY / SAVE
# ============================================================

def save_original_image(
    image,
    args,
):

    if not (
        args.display_original
        or args.output_origin
        or (args.save and args.display_original)
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
        save_requested_output(
            fig,
            args.output_origin,
        )

    if args.save and args.display_original:
        save_figure_png_pdf(
            fig,
            args,
            "original",
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

    if not (
        args.display_processed
        or args.output
        or (args.save and args.display_processed)
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

    # Background-subtracted image
    ax1.imshow(
        processed,
        origin="upper",
        norm=norm,
    )

    ax1.set_title(
        figure_title(
            args,
            "Background-subtracted data",
        )
    )

    ax1.axis("off")

    # Segmentation map
    ax2.imshow(
        segment_map.data,
        origin="upper",
        cmap=segment_map.cmap,
        interpolation="nearest",
    )

    ax2.set_title(
        figure_title(
            args,
            "Deblended Segmentation Image — All Detected Sources",
        )
    )

    ax2.axis("off")

    plt.tight_layout()

    if args.output:
        save_requested_output(
            fig,
            args.output,
        )

    if args.save and args.display_processed:
        save_figure_png_pdf(
            fig,
            args,
            "processed",
        )

    if args.display_processed:
        plt.show()

    plt.close(fig)


# ============================================================
# COORDINATE OUTPUT
# ============================================================

def save_python_coordinates(
    table,
    x_col,
    y_col,
    filename,
):
    """
    Save ALL detected-source coordinates for Python.

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
        " --- saving Python coordinates to:",
        filename,
    )

    with open(
        filename,
        "w",
        encoding="utf-8",
    ) as output_file:

        for row in table:

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
    table,
    x_col,
    y_col,
    image_height,
    filename,
):
    """
    Save ALL detected-source coordinates for ROOT.

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
        " --- saving ROOT coordinates to:",
        filename,
    )

    with open(
        filename,
        "w",
        encoding="utf-8",
    ) as output_file:

        for row in table:

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
# CIRCLE DISPLAY / SAVE
# ============================================================

def display_source_circles(
    processed,
    table,
    x_col,
    y_col,
    args,
):
    """
    Display ALL detected sources with the same fixed integration radius.
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
        figure_title(
            args,
            "Background subtracted data with detected sources"
            f"\nfixed integration radius = {integration_radius:.1f} px",
        )
    )

    for row in table:

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

    # Example:
    #   --circle --save
    # produces:
    #   <input_stem>_circles.png
    #   <input_stem>_circles.pdf
    save_figure_png_pdf(
        fig,
        args,
        "circles",
    )

    plt.show()
    plt.close(fig)


# ============================================================
# EMPTY OUTPUTS WHEN NO SOURCES ARE FOUND
# ============================================================

def create_empty_coordinate_files(args):

    if args.coordinates_python:

        Path(args.coordinates_python).parent.mkdir(
            parents=True,
            exist_ok=True,
        )

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

        Path(args.coordinates_root).parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        open(
            args.coordinates_root,
            "w",
            encoding="utf-8",
        ).close()

        print(
            " --- empty ROOT coordinate file created:",
            args.coordinates_root,
        )


# ============================================================
# MAIN PROGRAM
# ============================================================

def main():

    args = parse_arguments()

    # Validate configuration.
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
        " --- isolation selection = DISABLED "
        "(all detected/deblended sources are kept)"
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

    # Preserve the behavior of the NOisolated code:
    # convolution is used for detection only if --convolution is passed.
    if args.convolution:
        detection_image = convolved_data
        print(" --- detection image: Gaussian-convolved image")
    else:
        detection_image = processed
        print(" --- detection image: background-subtracted image")


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
        detection_image,
        threshold,
    )


    # ========================================================
    # No sources detected
    # ========================================================

    if segment_map is None:

        print(
            " --- No sources detected."
        )

        create_empty_coordinate_files(args)

        return 0


    # ========================================================
    # Source catalogue
    # ========================================================

    catalog = SourceCatalog(
        processed,
        segment_map,
        convolved_data=detection_image,
    )

    table = catalog.to_table()

    x_col, y_col = (
        get_centroid_column_names(table)
    )

    table[x_col].info.format = ".2f"
    table[y_col].info.format = ".2f"

    # One consecutive defect ID for EVERY detected/deblended source.
    table["defect_id"] = np.arange(
        1,
        len(table) + 1,
    )

    # Assign the SAME integration area/radius to every row.
    table["integration_area_pix2"] = np.full(
        len(table),
        fixed_area,
        dtype=float,
    )

    table["integration_radius_pix"] = np.full(
        len(table),
        integration_radius,
        dtype=float,
    )

    table["integration_area_pix2"].info.format = ".2f"
    table["integration_radius_pix"].info.format = ".2f"


    # ========================================================
    # Statistics / source list
    # ========================================================

    number_all = len(table)

    print(
        f" --- detected and deblended sources: "
        f"{number_all}"
    )

    print(
        " --- ALL detected sources are retained."
    )

    if number_all > 0:

        print(
            " --- source centroids:"
        )

        for row in table:

            print(
                f"     defect_id="
                f"{int(row['defect_id'])}, "
                f"x={to_float(row[x_col]):.2f}, "
                f"y={to_float(row[y_col]):.2f}, "
                f"R={integration_radius:.2f} px"
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
    # Focus on detected source
    # ========================================================

    if args.focus is not None:

        matches = np.where(
            np.array(
                table["defect_id"],
                dtype=int,
            )
            == args.focus
        )[0]

        if len(matches) == 0:

            print(
                f"ERROR: selected defect_id "
                f"{args.focus} is out of range. "
                f"Available values: 1 ... "
                f"{number_all}"
            )

            return 1


        focus_index = int(
            matches[0]
        )

        xc = to_float(
            table[x_col][focus_index]
        )

        yc = to_float(
            table[y_col][focus_index]
        )

        print(
            f" --- focus on defect "
            f"{args.focus}: "
            f"x={xc:.2f}, "
            f"y={yc:.2f}, "
            f"integration radius="
            f"{integration_radius:.2f} px"
        )

        display_focus(
            processed,
            xc,
            yc,
            integration_radius,
            (
                f"Defect {args.focus} "
                f"with radius "
                f"{integration_radius:.1f} px"
            ),
            args,
            f"focus_defect_{args.focus}",
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
                f"Manual focus: "
                f"x={xc:.1f}, "
                f"y={yc:.1f}, "
                f"r={manual_radius:.1f} px"
            ),
            args,
            (
                f"focus_area_"
                f"x{xc:.1f}_y{yc:.1f}_r{manual_radius:.1f}"
                .replace(".", "p")
            ),
        )


    # ========================================================
    # Save coordinates
    # ========================================================

    if args.coordinates_python:

        Path(args.coordinates_python).parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        save_python_coordinates(
            table,
            x_col,
            y_col,
            args.coordinates_python,
        )


    if args.coordinates_root:

        Path(args.coordinates_root).parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        save_root_coordinates(
            table,
            x_col,
            y_col,
            processed.shape[0],
            args.coordinates_root,
        )


    # ========================================================
    # Display / save circles
    # ========================================================

    if args.circles:

        display_source_circles(
            processed,
            table,
            x_col,
            y_col,
            args,
        )


    return 0


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    sys.exit(main())
