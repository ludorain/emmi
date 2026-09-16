#!/usr/bin/env python3

"""
Create a 2x2 validation display for one global hotspot:
  top-left:  before_annealing reference image (already rotated)
  top-right: target annealing image before rotation/shift (processed)
  bottom-left: same reference image
  bottom-right: target annealing image after rotation/shift

All panels show a 50x50 px crop centred on the reference global coordinates.
The same small red X is drawn in all panels, so misalignment is visible before
rotation/shift and corrected after alignment.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.background import Background2D, MedianBackground


# Exact display configuration used by the original --focus implementation.
vmin_global = 0.0
vmax_global = 5.0
background_box_size = (50, 50)
background_filter_size = (3, 3)


def read_image(path: Path) -> np.ndarray:
    image = np.asarray(tifffile.imread(path), dtype=float)
    if image.ndim != 2:
        raise ValueError(f"Expected a 2D TIF image, got shape {image.shape} for {path}")
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


def crop_focus(image: np.ndarray, xc: float, yc: float, crop_size: int):
    """
    Crop a fixed-size square around the selected coordinates, padding with NaN
    when the requested area exceeds the image boundaries.
    """
    half = crop_size // 2
    x0 = int(round(xc)) - half
    y0 = int(round(yc)) - half
    x1 = x0 + crop_size
    y1 = y0 + crop_size

    crop = np.full((crop_size, crop_size), np.nan, dtype=float)

    src_x0 = max(0, x0)
    src_y0 = max(0, y0)
    src_x1 = min(image.shape[1], x1)
    src_y1 = min(image.shape[0], y1)

    if src_x0 < src_x1 and src_y0 < src_y1:
        dst_x0 = src_x0 - x0
        dst_y0 = src_y0 - y0
        dst_x1 = dst_x0 + (src_x1 - src_x0)
        dst_y1 = dst_y0 + (src_y1 - src_y0)
        crop[dst_y0:dst_y1, dst_x0:dst_x1] = image[src_y0:src_y1, src_x0:src_x1]

    x_local = xc - x0
    y_local = yc - y0
    return crop, x_local, y_local


def focus_norm():
    """Exact normalization used by the original --focus display."""
    return ImageNormalize(
        vmin=vmin_global,
        vmax=vmax_global,
        stretch=SqrtStretch(),
    )


def draw_x(ax, x: float, y: float, half_arm: float = 2.8, color: str = "red"):
    ax.plot([x - half_arm, x + half_arm], [y - half_arm, y + half_arm],
            color=color, linewidth=1.6, solid_capstyle="round")
    ax.plot([x - half_arm, x + half_arm], [y + half_arm, y - half_arm],
            color=color, linewidth=1.6, solid_capstyle="round")

def draw_crosshair(ax, x: float, y: float, radius: float = 5.0, arm_len: float = 7.5, color: str = "white"):
    """Draw a crosshair centered at (x, y) with a circle and cross lines."""
    # Cerchio centrale
    circle = plt.Circle((x, y), radius=radius, color=color, fill=False, linewidth=1.2)
    ax.add_patch(circle)

    # Linea orizzontale
    ax.plot([x - arm_len, x + arm_len], [y, y], color=color, linewidth=1.2)
    # Linea verticale
    ax.plot([x, x], [y - arm_len, y + arm_len], color=color, linewidth=1.2)

def draw_panel(ax, crop, x_local, y_local, title, norm):
    ax.imshow(crop, origin="upper", norm=norm)
    draw_crosshair(ax, x_local, y_local)
    ax.set_title(title, fontsize=10)
    ax.set_xlim(0, crop.shape[1])
    ax.set_ylim(crop.shape[0], 0)
    ax.set_aspect("equal")
    ax.axis("off")


def make_four_panel_figure(reference_image, target_before, target_after,
                           xc, yc, crop_size, suptitle, output_path):
    ref_crop, x_local, y_local = crop_focus(reference_image, xc, yc, crop_size)
    before_crop, _, _ = crop_focus(target_before, xc, yc, crop_size)
    after_crop, _, _ = crop_focus(target_after, xc, yc, crop_size)

    crops = [ref_crop, before_crop, ref_crop, after_crop]
    norm = focus_norm()

    fig, axes = plt.subplots(2, 2, figsize=(8.4, 8.2))

    draw_panel(axes[0, 0], ref_crop, x_local, y_local,
               "Reference: before annealing (rotated)", norm)
    draw_panel(axes[0, 1], before_crop, x_local, y_local,
               "Selected phase: before rotation + shift", norm)
    draw_panel(axes[1, 0], ref_crop, x_local, y_local,
               "Reference: before annealing (rotated)", norm)
    draw_panel(axes[1, 1], after_crop, x_local, y_local,
               "Selected phase: after rotation + shift", norm)

    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Create a 4-panel rotation/shift validation display for one global hotspot."
    )
    parser.add_argument("--reference", required=True, type=Path,
                        help="before_annealing data=diff image in 3rotated")
    parser.add_argument("--target-before", required=True, type=Path,
                        help="selected-phase data=diff image in 2processed")
    parser.add_argument("--target-after", required=True, type=Path,
                        help="selected-phase data=diff image in 3rotated")
    parser.add_argument("--x", required=True, type=float,
                        help="reference hotspot x coordinate in ROOT convention")
    parser.add_argument("--y", required=True, type=float,
                        help="reference hotspot y coordinate in ROOT convention")
    parser.add_argument("--sensor", required=True)
    parser.add_argument("--constant-label", required=True,
                        help="e.g. T=20 or v=5")
    parser.add_argument("--phase", required=True)
    parser.add_argument("--spot", required=True, type=int)
    parser.add_argument("--scan-label", required=True,
                        help="operating point used for the displayed data=diff image")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--crop-size", type=int, default=50)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    if args.crop_size <= 0:
        raise SystemExit("--crop-size must be positive")

    for path in (args.reference, args.target_before, args.target_after):
        if not path.is_file():
            raise FileNotFoundError(path)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    reference = background_subtract(read_image(args.reference))
    target_before = background_subtract(read_image(args.target_before))
    target_after = background_subtract(read_image(args.target_after))

    # Convert coordinates from ROOT convention to Python/image convention
    # ROOT file stores:
    #   x_root = x_python
    #   y_root = image_height - y_python
    # therefore:
    #   y_python = image_height - y_root
    x_python = args.x
    y_python = reference.shape[0] - args.y

    phase_safe = args.phase.replace("/", "_")
    common = f"{args.sensor}_{args.constant_label}_{phase_safe}_spot{args.spot}"
    output_path = args.output_dir / f"{common}_rotation_shift_validation.png"

    suptitle = (
        f"{args.sensor} | global hotspot {args.spot} | {args.constant_label} | "
        f"{args.phase} | {args.scan_label} | reference center=({x_python:.2f}, {y_python:.2f})"
    )

    make_four_panel_figure(
        reference, target_before, target_after,
        x_python, y_python, args.crop_size,
        suptitle, output_path,
    )

    print(f"Created: {output_path}")

    if args.show:
        image = plt.imread(output_path)
        plt.figure(figsize=(9, 9))
        plt.imshow(image)
        plt.axis("off")
        plt.show()


if __name__ == "__main__":
    main()
