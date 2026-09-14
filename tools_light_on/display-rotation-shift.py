#!/usr/bin/env python3

"""
Display the same global EMMI hotspot before and after geometrical alignment.

The program receives three already-selected TIF images:
  1) before_annealing reference image after rotation;
  2) target annealing image after cleanup but before rotation/shift;
  3) target annealing image after the final rotation/shift.

A fixed 60x60 px crop, centred on the reference global coordinates, is shown
for all panels. The same reference crosshair is drawn in every image.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile
from matplotlib.patches import Circle


def read_image(path: Path) -> np.ndarray:
    image = np.asarray(tifffile.imread(path), dtype=float)
    if image.ndim != 2:
        raise ValueError(f"Expected a 2D TIF image, got shape {image.shape} for {path}")
    return image


def crop_fixed_size(image: np.ndarray, x: float, y: float, size: int):
    """Return an exactly size x size crop centred on (x, y), padding if needed."""
    half = size // 2
    cx = int(round(x))
    cy = int(round(y))

    x0 = cx - half
    y0 = cy - half
    x1 = x0 + size
    y1 = y0 + size

    crop = np.full((size, size), np.nan, dtype=float)

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

    # Exact sub-pixel reference position inside the crop.
    cross_x = x - x0
    cross_y = y - y0
    return crop, cross_x, cross_y


def robust_limits(crop: np.ndarray):
    finite = crop[np.isfinite(crop)]
    if finite.size == 0:
        return 0.0, 1.0

    vmin, vmax = np.percentile(finite, [1.0, 99.7])
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin = float(np.min(finite))
        vmax = float(np.max(finite))
    if vmax <= vmin:
        vmax = vmin + 1.0
    return vmin, vmax


def draw_crosshair(ax, x: float, y: float):
    """Draw a small red target without covering the hotspot centre."""
    inner = 4.0
    outer = 11.0

    ax.plot([x - outer, x - inner], [y, y], color="red", linewidth=1.4)
    ax.plot([x + inner, x + outer], [y, y], color="red", linewidth=1.4)
    ax.plot([x, x], [y - outer, y - inner], color="red", linewidth=1.4)
    ax.plot([x, x], [y + inner, y + outer], color="red", linewidth=1.4)
    ax.add_patch(Circle((x, y), radius=3.2, fill=False, edgecolor="red", linewidth=1.2))


def draw_panel(ax, image, x, y, crop_size, title):
    crop, cross_x, cross_y = crop_fixed_size(image, x, y, crop_size)
    vmin, vmax = robust_limits(crop)

    ax.imshow(crop, cmap="gray", origin="upper", vmin=vmin, vmax=vmax)
    draw_crosshair(ax, cross_x, cross_y)
    ax.set_title(title, fontsize=11)
    ax.set_xlim(-0.5, crop_size - 0.5)
    ax.set_ylim(crop_size - 0.5, -0.5)
    ax.set_aspect("equal")
    ax.axis("off")


def make_pair_figure(reference_image, target_image, x, y, crop_size,
                     left_title, right_title, suptitle, output_path):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.3))

    draw_panel(axes[0], reference_image, x, y, crop_size, left_title)
    draw_panel(axes[1], target_image, x, y, crop_size, right_title)

    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Compare an EMMI hotspot before and after rotation/shift alignment."
    )
    parser.add_argument("--reference", required=True, type=Path,
                        help="before_annealing data=diff image in 3rotated")
    parser.add_argument("--target-before", required=True, type=Path,
                        help="target-phase data=diff image in 2processed")
    parser.add_argument("--target-after", required=True, type=Path,
                        help="target-phase data=diff image in 3rotated")
    parser.add_argument("--x", required=True, type=float,
                        help="reference hotspot x coordinate")
    parser.add_argument("--y", required=True, type=float,
                        help="reference hotspot y coordinate")
    parser.add_argument("--sensor", required=True)
    parser.add_argument("--constant-label", required=True,
                        help="e.g. T=20 or v=5")
    parser.add_argument("--phase", required=True)
    parser.add_argument("--spot", required=True, type=int)
    parser.add_argument("--scan-label", required=True,
                        help="operating point used for the displayed data=diff image")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--crop-size", type=int, default=60)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    if args.crop_size <= 0:
        raise SystemExit("--crop-size must be positive")

    for path in (args.reference, args.target_before, args.target_after):
        if not path.is_file():
            raise FileNotFoundError(path)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    reference = read_image(args.reference)
    target_before = read_image(args.target_before)
    target_after = read_image(args.target_after)

    phase_safe = args.phase.replace("/", "_")
    common = f"{args.sensor}_{args.constant_label}_{phase_safe}_spot{args.spot}"

    before_output = args.output_dir / f"{common}_before_rotation_shift.png"
    after_output = args.output_dir / f"{common}_after_rotation_shift.png"

    suptitle = (
        f"{args.sensor} | global hotspot {args.spot} | {args.constant_label} | "
        f"{args.phase} | {args.scan_label} | reference center=({args.x:.2f}, {args.y:.2f})"
    )

    make_pair_figure(
        reference,
        target_before,
        args.x,
        args.y,
        args.crop_size,
        "Reference: before annealing (rotated)",
        f"{args.phase}: before rotation + shift",
        suptitle,
        before_output,
    )

    make_pair_figure(
        reference,
        target_after,
        args.x,
        args.y,
        args.crop_size,
        "Reference: before annealing (rotated)",
        f"{args.phase}: after rotation + shift",
        suptitle,
        after_output,
    )

    print(f"Created: {before_output}")
    print(f"Created: {after_output}")

    if args.show:
        # Re-open the saved images only when explicitly requested.
        for output in (before_output, after_output):
            image = plt.imread(output)
            plt.figure(figsize=(10, 5))
            plt.imshow(image)
            plt.axis("off")
            plt.show()


if __name__ == "__main__":
    main()
