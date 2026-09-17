#!/usr/bin/env python3

"""
Create a 1x3 validation display for one global hotspot:

  left:   before_annealing reference image
  center: target annealing image before rotation/shift
  right:  target annealing image after rotation/shift

All panels show the same crop centred on the reference global coordinates.
The same crosshair is drawn in all panels, so misalignment is visible before
rotation/shift and corrected after alignment.

The figure is saved both as a high-resolution PNG and as a PDF.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.background import Background2D, MedianBackground


# ============================================================
# DISPLAY CONFIGURATION
# ============================================================

vmin_global = 0.0
vmax_global = 5.0

background_box_size = (50, 50)
background_filter_size = (3, 3)


# ============================================================
# IMAGE READING / PROCESSING
# ============================================================

def read_image(path: Path) -> np.ndarray:
    image = np.asarray(tifffile.imread(path), dtype=float)

    if image.ndim != 2:
        raise ValueError(
            f"Expected a 2D TIF image, got shape {image.shape} for {path}"
        )

    return image


def background_subtract(image: np.ndarray) -> np.ndarray:
    """Replicate the processing used before display_focus in the source code."""

    background = Background2D(
        image,
        box_size=background_box_size,
        filter_size=background_filter_size,
        bkg_estimator=MedianBackground(),
    )

    return image - background.background


# ============================================================
# CROP
# ============================================================

def crop_focus(
    image: np.ndarray,
    xc: float,
    yc: float,
    crop_size: int
):
    """
    Crop a fixed-size square around the selected coordinates.

    If the requested crop extends outside the image boundaries,
    the missing region is padded with NaN.
    """

    half = crop_size // 2

    x0 = int(round(xc)) - half
    y0 = int(round(yc)) - half

    x1 = x0 + crop_size
    y1 = y0 + crop_size

    crop = np.full(
        (crop_size, crop_size),
        np.nan,
        dtype=float
    )

    src_x0 = max(0, x0)
    src_y0 = max(0, y0)

    src_x1 = min(image.shape[1], x1)
    src_y1 = min(image.shape[0], y1)

    if src_x0 < src_x1 and src_y0 < src_y1:

        dst_x0 = src_x0 - x0
        dst_y0 = src_y0 - y0

        dst_x1 = dst_x0 + (src_x1 - src_x0)
        dst_y1 = dst_y0 + (src_y1 - src_y0)

        crop[
            dst_y0:dst_y1,
            dst_x0:dst_x1
        ] = image[
            src_y0:src_y1,
            src_x0:src_x1
        ]

    x_local = xc - x0
    y_local = yc - y0

    return crop, x_local, y_local


# ============================================================
# DISPLAY
# ============================================================

def focus_norm():
    """Normalization used by the original --focus display."""

    return ImageNormalize(
        vmin=vmin_global,
        vmax=vmax_global,
        stretch=SqrtStretch(),
    )


def draw_crosshair(
    ax,
    x: float,
    y: float,
    radius: float = 5.0,
    arm_len: float = 7.5,
    color: str = "white",
):
    """Draw a crosshair centered at (x, y)."""

    circle = plt.Circle(
        (x, y),
        radius=radius,
        color=color,
        fill=False,
        linewidth=1.2,
    )

    ax.add_patch(circle)

    # Horizontal line
    ax.plot(
        [x - arm_len, x + arm_len],
        [y, y],
        color=color,
        linewidth=1.2,
    )

    # Vertical line
    ax.plot(
        [x, x],
        [y - arm_len, y + arm_len],
        color=color,
        linewidth=1.2,
    )


def draw_panel(
    ax,
    crop,
    x_local,
    y_local,
    title,
    norm,
):
    ax.imshow(
        crop,
        origin="upper",
        norm=norm,
        interpolation="none",
    )

    draw_crosshair(
        ax,
        x_local,
        y_local
    )

    ax.set_title(
        title,
        fontsize=11
    )

    ax.set_xlim(
        0,
        crop.shape[1]
    )

    ax.set_ylim(
        crop.shape[0],
        0
    )

    ax.set_aspect("equal")
    ax.axis("off")


# ============================================================
# FIGURE
# ============================================================

def make_three_panel_figure(
    reference_image,
    target_before,
    target_after,
    xc,
    yc,
    crop_size,
    suptitle,
    output_base,
):

    # --------------------------------------------------------
    # Create crops
    # --------------------------------------------------------

    ref_crop, x_local, y_local = crop_focus(
        reference_image,
        xc,
        yc,
        crop_size
    )

    before_crop, _, _ = crop_focus(
        target_before,
        xc,
        yc,
        crop_size
    )

    after_crop, _, _ = crop_focus(
        target_after,
        xc,
        yc,
        crop_size
    )

    norm = focus_norm()

    # --------------------------------------------------------
    # 1 x 3 layout
    # --------------------------------------------------------

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(12.6, 4.5)
    )

    draw_panel(
        axes[0],
        ref_crop,
        x_local,
        y_local,
        "Reference: before annealing",
        norm,
    )

    draw_panel(
        axes[1],
        before_crop,
        x_local,
        y_local,
        "Selected phase: before rotation + shift",
        norm,
    )

    draw_panel(
        axes[2],
        after_crop,
        x_local,
        y_local,
        "Selected phase: after rotation + shift",
        norm,
    )

    fig.suptitle(
        suptitle,
        fontsize=12
    )

    fig.tight_layout(
        rect=(0, 0, 1, 0.92)
    )

    # --------------------------------------------------------
    # Save PNG
    # --------------------------------------------------------

    png_path = output_base.with_suffix(".png")

    fig.savefig(
        png_path,
        dpi=600,
        bbox_inches="tight",
    )

    # --------------------------------------------------------
    # Save PDF
    # --------------------------------------------------------

    pdf_path = output_base.with_suffix(".pdf")

    fig.savefig(
        pdf_path,
        format="pdf",
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(fig)

    return png_path, pdf_path


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Create a 3-panel rotation/shift validation display "
            "for one global hotspot."
        )
    )

    parser.add_argument(
        "--reference",
        required=True,
        type=Path,
        help="before_annealing data=diff image in 3rotated",
    )

    parser.add_argument(
        "--target-before",
        required=True,
        type=Path,
        help="selected-phase data=diff image in 2processed",
    )

    parser.add_argument(
        "--target-after",
        required=True,
        type=Path,
        help="selected-phase data=diff image in 3rotated",
    )

    parser.add_argument(
        "--x",
        required=True,
        type=float,
        help="reference hotspot x coordinate in ROOT convention",
    )

    parser.add_argument(
        "--y",
        required=True,
        type=float,
        help="reference hotspot y coordinate in ROOT convention",
    )

    parser.add_argument(
        "--sensor",
        required=True
    )

    parser.add_argument(
        "--constant-label",
        required=True,
        help="e.g. T=20 or v=5",
    )

    parser.add_argument(
        "--phase",
        required=True
    )

    parser.add_argument(
        "--spot",
        required=True,
        type=int
    )

    parser.add_argument(
        "--scan-label",
        required=True,
        help="operating point used for the displayed data=diff image",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path
    )

    parser.add_argument(
        "--crop-size",
        type=int,
        default=50
    )

    parser.add_argument(
        "--show",
        action="store_true"
    )

    args = parser.parse_args()

    if args.crop_size <= 0:
        raise SystemExit(
            "--crop-size must be positive"
        )

    for path in (
        args.reference,
        args.target_before,
        args.target_after,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Read and background-subtract images
    # --------------------------------------------------------

    reference = background_subtract(
        read_image(args.reference)
    )

    target_before = background_subtract(
        read_image(args.target_before)
    )

    target_after = background_subtract(
        read_image(args.target_after)
    )

    # --------------------------------------------------------
    # Convert ROOT coordinates to image/Python convention
    #
    # x_root = x_python
    # y_root = image_height - y_python
    #
    # therefore:
    #
    # y_python = image_height - y_root
    # --------------------------------------------------------

    x_python = args.x
    y_python = reference.shape[0] - args.y

    # --------------------------------------------------------
    # Output filename
    # --------------------------------------------------------

    phase_safe = args.phase.replace("/", "_")

    common = (
        f"{args.sensor}_"
        f"{args.constant_label}_"
        f"{phase_safe}_"
        f"spot{args.spot}_"
        f"rotation_shift_validation"
    )

    output_base = (
        args.output_dir / common
    )

    # --------------------------------------------------------
    # Figure title
    # --------------------------------------------------------

    suptitle = (
        f"{args.sensor} | global hotspot {args.spot} | "
        f"{args.constant_label} | "
        f"{args.phase} | "
        f"{args.scan_label} | "
        f"reference center=({x_python:.2f}, {y_python:.2f})"
    )

    # --------------------------------------------------------
    # Create figure
    # --------------------------------------------------------

    png_path, pdf_path = make_three_panel_figure(
        reference,
        target_before,
        target_after,
        x_python,
        y_python,
        args.crop_size,
        suptitle,
        output_base,
    )

    print(f"Created PNG: {png_path}")
    print(f"Created PDF: {pdf_path}")

    # --------------------------------------------------------
    # Optional display
    # --------------------------------------------------------

    if args.show:

        image = plt.imread(png_path)

        plt.figure(
            figsize=(13, 5)
        )

        plt.imshow(image)

        plt.axis("off")
        plt.show()


if __name__ == "__main__":
    main()