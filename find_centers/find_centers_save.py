#!/usr/bin/env python3
"""
Detect isolated EMMI hotspots and save a focused crop.

The source-finding logic is intentionally kept identical to
find_defects_isolated_changeR.py:
  - Background2D with MedianBackground
  - background box (150, 150), filter (3, 3)
  - threshold = 3 * background RMS
  - Gaussian kernel FWHM=3 px, size=5
  - SourceFinder(npixels=30, deblend=True, nlevels=32, contrast=0.001)
  - variable equivalent radius from the segmented area
  - isolation cut: >60 px from every other source and from every border

Two focus modes are available:
  --focus LOCAL_SPOT
      LOCAL_SPOT is the 0-based local hotspot ID used in the mapping CSV.
      The crop is centred on the centroid found in this image.

  --focus_area X Y
      X,Y are Python image coordinates. They are used directly as the crop
      centre. This mode is intended for phases in which the global hotspot is
      not detected; the pipeline passes the coordinates of its last detection.

In both modes the red circle always has radius 20 pixels.
"""

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile

from astropy.convolution import convolve
from astropy.visualization import ImageNormalize, SqrtStretch
from matplotlib.patches import Circle
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import SourceCatalog, SourceFinder, make_2dgaussian_kernel


# ============================================================
# SOURCE-FINDING CONFIGURATION
# Identical to find_defects_isolated_changeR.py
# ============================================================
isolation_radius = 60.0
background_box_size = (150, 150)
background_filter_size = (3, 3)
threshold_sigma = 3.0
kernel_fwhm = 3.0
kernel_size = 5

# ============================================================
# FOCUS DISPLAY CONFIGURATION
# ============================================================
CIRCLE_RADIUS = 20.0
CROP_HALF_SIZE = 50  # same fixed half-size used by find_centers_focus_area.py


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Detect isolated luminous hotspots and save a focused crop with "
            "a fixed 20-pixel circle."
        )
    )

    parser.add_argument("--input", required=True, help="Input 2D TIF image")

    focus_group = parser.add_mutually_exclusive_group(required=True)
    focus_group.add_argument(
        "--focus",
        type=int,
        help=(
            "0-based local hotspot ID, matching the local_spot column of the "
            "global mapping CSV"
        ),
    )
    focus_group.add_argument(
        "--focus_area",
        nargs=2,
        type=float,
        metavar=("X", "Y"),
        help="Manual focus position in Python image coordinates",
    )

    parser.add_argument(
        "--output_focus",
        default=None,
        help="Output basename; .png and .pdf are appended",
    )
    parser.add_argument(
        "--focus_title",
        default="Focused hotspot",
        help="Title written above the focused image",
    )

    # Shared display scale supplied by hotspot_photos.sh.
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument(
        "--measure_max",
        action="store_true",
        help=(
            "Print the maximum pixel value in the focused background-subtracted "
            "crop. Used by the pipeline to calculate one common scale per scan."
        ),
    )

    # Keep SourceFinder options from the original code.
    parser.add_argument("--npixels", type=int, default=30)
    parser.add_argument("--nlevels", type=int, default=32)
    parser.add_argument("--contrast", type=float, default=0.001)

    return parser.parse_args()


def validate_configuration(args):
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
    if args.focus is not None and args.focus < 0:
        raise ValueError("--focus/local_spot must be >= 0")


def to_float(value):
    if hasattr(value, "value"):
        value = value.value
    return float(value)


def get_centroid_column_names(table):
    if "x_centroid" in table.colnames and "y_centroid" in table.colnames:
        return "x_centroid", "y_centroid"
    if "xcentroid" in table.colnames and "ycentroid" in table.colnames:
        return "xcentroid", "ycentroid"
    raise RuntimeError("Centroid columns were not found in the SourceCatalog table.")


def calculate_source_dimensions(catalog, table):
    segment_area = np.array(
        [to_float(area) for area in catalog.segment_area], dtype=float
    )
    equivalent_radius = np.sqrt(segment_area / np.pi)
    table["segment_area_pix2"] = segment_area
    table["equivalent_radius_pix"] = equivalent_radius
    table["segment_area_pix2"].info.format = ".2f"
    table["equivalent_radius_pix"].info.format = ".2f"


def select_isolated_sources(table, x_col, y_col, image_shape):
    """Copied in logic from find_defects_isolated_changeR.py."""
    ny, nx = image_shape

    x = np.array([to_float(value) for value in table[x_col]], dtype=float)
    y = np.array([to_float(value) for value in table[y_col]], dtype=float)
    number_of_sources = len(table)

    if number_of_sources == 1:
        nearest_source_distance = np.array([np.inf], dtype=float)
    else:
        coordinates = np.column_stack((x, y))
        differences = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
        distance_matrix = np.sqrt(np.sum(differences**2, axis=2))
        np.fill_diagonal(distance_matrix, np.inf)
        nearest_source_distance = np.min(distance_matrix, axis=1)

    isolated_from_sources = nearest_source_distance > isolation_radius

    distance_left = x
    distance_right = (nx - 1) - x
    distance_top = y
    distance_bottom = (ny - 1) - y

    border_distance = np.minimum.reduce(
        [distance_left, distance_right, distance_top, distance_bottom]
    )
    isolated_from_border = border_distance > isolation_radius

    isolated_mask = isolated_from_sources & isolated_from_border

    isolated_table = table[isolated_mask].copy()
    isolated_table["nearest_source_distance_pix"] = nearest_source_distance[
        isolated_mask
    ]
    isolated_table["nearest_border_distance_pix"] = border_distance[isolated_mask]

    # IMPORTANT: the mapping pipeline uses 0-based local_spot IDs.
    isolated_table["local_spot"] = np.arange(len(isolated_table), dtype=int)

    return (
        isolated_table,
        isolated_mask,
        nearest_source_distance,
        border_distance,
    )


def read_and_process_image(filename, args):
    image = tifffile.imread(filename).astype(float)
    if image.ndim != 2:
        raise ValueError(f"expected a 2D image, got shape {image.shape}")

    background_estimator = MedianBackground()
    background = Background2D(
        image,
        box_size=background_box_size,
        filter_size=background_filter_size,
        bkg_estimator=background_estimator,
    )

    processed = image - background.background
    threshold = threshold_sigma * background.background_rms

    kernel = make_2dgaussian_kernel(kernel_fwhm, size=kernel_size)
    convolved_data = convolve(processed, kernel)

    finder = SourceFinder(
        npixels=args.npixels,
        deblend=True,
        nlevels=args.nlevels,
        contrast=args.contrast,
        progress_bar=False,
    )

    segment_map = finder(convolved_data, threshold)

    if segment_map is None:
        return image, processed, None, None, None, None

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

    (
        isolated_table,
        _,
        _,
        _,
    ) = select_isolated_sources(source_table, x_col, y_col, processed.shape)

    return image, processed, isolated_table, x_col, y_col, segment_map


def resolve_focus_position(args, processed, isolated_table, x_col, y_col):
    if args.focus is not None:
        if isolated_table is None or len(isolated_table) == 0:
            raise ValueError(
                "--focus was requested, but no isolated hotspots were detected"
            )

        matches = np.where(
            np.asarray(isolated_table["local_spot"], dtype=int) == args.focus
        )[0]

        if len(matches) != 1:
            available = [int(v) for v in isolated_table["local_spot"]]
            raise ValueError(
                f"local hotspot {args.focus} was not found; available local_spot "
                f"IDs are {available}"
            )

        focus_index = int(matches[0])
        xc = to_float(isolated_table[x_col][focus_index])
        yc = to_float(isolated_table[y_col][focus_index])

        print(
            f" --- focus local_spot={args.focus}: "
            f"centroid x={xc:.3f}, y={yc:.3f} [Python coordinates]"
        )
        return xc, yc

    # --focus_area: coordinates are already Python coordinates.
    xc, yc = args.focus_area
    print(
        f" --- focus_area: x={xc:.3f}, y={yc:.3f} [Python coordinates]"
    )
    return float(xc), float(yc)


def make_focus_crop(processed, xc, yc):
    """Same crop-centering logic as find_centers_focus_area.py."""
    ny, nx = processed.shape

    if not (math.isfinite(xc) and math.isfinite(yc)):
        raise ValueError(f"non-finite focus coordinate x={xc}, y={yc}")
    if not (0.0 <= xc < nx and 0.0 <= yc < ny):
        raise ValueError(
            f"focus coordinate x={xc:.3f}, y={yc:.3f} is outside image "
            f"size {nx}x{ny}"
        )

    x_min = max(0, int(round(xc)) - CROP_HALF_SIZE)
    x_max = min(nx, int(round(xc)) + CROP_HALF_SIZE)
    y_min = max(0, int(round(yc)) - CROP_HALF_SIZE)
    y_max = min(ny, int(round(yc)) + CROP_HALF_SIZE)

    if x_min >= x_max or y_min >= y_max:
        raise ValueError(
            f"empty crop for x={xc:.3f}, y={yc:.3f}; check coordinates"
        )

    focus_image = processed[y_min:y_max, x_min:x_max]
    xc_local = xc - x_min
    yc_local = yc - y_min

    return focus_image, xc_local, yc_local


def save_focus_image(focus_image, xc_local, yc_local, args):
    if args.output_focus is None:
        raise ValueError("--output_focus is required unless --measure_max is used")
    if args.vmin is None or args.vmax is None:
        raise ValueError("--vmin and --vmax are required when saving")
    if not (
        math.isfinite(args.vmin)
        and math.isfinite(args.vmax)
        and args.vmax > args.vmin
    ):
        raise ValueError(f"invalid display limits: {args.vmin}, {args.vmax}")

    norm = ImageNormalize(
        vmin=args.vmin,
        vmax=args.vmax,
        stretch=SqrtStretch(),
    )

    fig_focus, ax_focus = plt.subplots(figsize=(7.4, 7.6))
    ax_focus.imshow(focus_image, origin="upper", norm=norm)

    # Always draw a fixed R=20 pixel circle at the exact focus coordinate.
    circle = Circle(
        (xc_local, yc_local),
        radius=CIRCLE_RADIUS,
        edgecolor="red",
        facecolor="none",
        linewidth=1.5,
    )
    ax_focus.add_patch(circle)

    ax_focus.set_title(args.focus_title, fontsize=14, pad=14, wrap=True)
    ax_focus.set_xlim(0, focus_image.shape[1])
    ax_focus.set_ylim(focus_image.shape[0], 0)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    output_base = Path(args.output_focus)
    output_base.parent.mkdir(parents=True, exist_ok=True)

    png_name = str(output_base) + ".png"
    pdf_name = str(output_base) + ".pdf"

    fig_focus.savefig(png_name, dpi=300, bbox_inches="tight")
    fig_focus.savefig(pdf_name, bbox_inches="tight")
    plt.close(fig_focus)

    print(
        f" --- red circle centre in crop: x={xc_local:.3f}, y={yc_local:.3f}; "
        f"R={CIRCLE_RADIUS:.1f} px"
    )
    print(f" --- saved: {png_name}")
    print(f" --- saved: {pdf_name}")


def main():
    args = parse_arguments()

    try:
        validate_configuration(args)

        print(" --- opening input image:", args.input)
        (
            image,
            processed,
            isolated_table,
            x_col,
            y_col,
            segment_map,
        ) = read_and_process_image(args.input, args)

        if isolated_table is None:
            print(" --- no sources detected by SourceFinder")
        else:
            print(f" --- isolated hotspots found: {len(isolated_table)}")

        xc, yc = resolve_focus_position(
            args, processed, isolated_table, x_col, y_col
        )
        focus_image, xc_local, yc_local = make_focus_crop(processed, xc, yc)

        if args.measure_max:
            finite = focus_image[np.isfinite(focus_image)]
            if finite.size == 0:
                raise ValueError("focused crop contains no finite pixels")
            # Final stdout line intentionally numeric for the shell pipeline.
            print(f"{float(np.max(finite)):.12g}")
            return 0

        save_focus_image(focus_image, xc_local, yc_local, args)
        return 0

    except (
        OSError,
        ValueError,
        RuntimeError,
        tifffile.TiffFileError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
