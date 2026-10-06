#!/usr/bin/env python3
"""
Detect luminous hotspots exactly as in find_defects_NOisolated_changeR.py and
save a focused crop.

Main features:
- source finding logic copied from find_defects_NOisolated_changeR.py
- no isolation cut
- --focus DEFECT_ID uses the detected hotspot list (1-based IDs as in the
  original script)
- --focus_area X Y RADIUS focuses on a manual area, following the logic of
  find_centers_focus_area.py
- the red circle is always drawn and has radius 20 pixels when saving the crop
- optional --detected_csv writes the detected hotspots for pipeline matching
- optional --measure_max prints the maximum value inside the focused crop, so a
  scan can share one common colour scale
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

CIRCLE_RADIUS = 20.0
FOCUS_HALF_SIZE = 40


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process an EMMI image')

    parser.add_argument('--input', type=str, required=True, help='Input TIF filename')

    focus_group = parser.add_mutually_exclusive_group(required=False)
    focus_group.add_argument('--focus', type=int, required=False,
                             help='Show a region around the selected defect_id written in --detected_csv (1-based)')
    focus_group.add_argument('--focus_area', nargs=3, type=float, metavar=('X', 'Y', 'RADIUS'),
                             help='Show a region around manually specified image coordinates: --focus_area x y radius')

    parser.add_argument('--output_focus', type=str, required=False,
                        help='Output basename for the focused image; both PNG and PDF are saved')
    parser.add_argument('--focus_title', type=str, default='Focused area',
                        help='Title of the focused figure')

    parser.add_argument('--detected_csv', type=str, required=False,
                        help='Optional output CSV for detected hotspots (defect_id,x,y,segment_area_pix2,equivalent_radius_pix)')

    parser.add_argument('--measure_max', action='store_true',
                        help='Measure and print the maximum value in the focused crop')
    parser.add_argument('--vmin', type=float, required=False,
                        help='Optional display vmin for the focus image')
    parser.add_argument('--vmax', type=float, required=False,
                        help='Optional display vmax for the focus image')

    parser.add_argument('--display_original', action='store_true', help='Display original image')
    parser.add_argument('--display_processed', action='store_true', help='Display processed image')
    parser.add_argument('--convolution', action='store_true', help='Apply convolution to the image for source detection')

    parser.add_argument('--npixels', type=int, default=36,
                        help='Minimum number of connected pixels for source detection')
    parser.add_argument('--nlevels', type=int, default=32,
                        help='Number of multi-thresholding levels for deblending')
    parser.add_argument('--contrast', type=float, default=0.001,
                        help='Deblending contrast parameter')

    return parser.parse_args()


def get_centroid_column_names(tbl):
    if 'x_centroid' in tbl.colnames and 'y_centroid' in tbl.colnames:
        return 'x_centroid', 'y_centroid'
    if 'xcentroid' in tbl.colnames and 'ycentroid' in tbl.colnames:
        return 'xcentroid', 'ycentroid'
    raise RuntimeError('Centroid columns not found in SourceCatalog table.')


def to_float(value):
    if hasattr(value, 'value'):
        value = value.value
    return float(value)


def save_detected_csv(tbl, x_col, y_col, filename):
    output_path = Path(filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w', encoding='utf-8') as f:
        f.write('defect_id,x,y,segment_area_pix2,equivalent_radius_pix\n')
        for row in tbl:
            defect_id = int(row['defect_id'])
            x = to_float(row[x_col])
            y = to_float(row[y_col])
            area = to_float(row['segment_area_pix2'])
            radius = to_float(row['equivalent_radius_pix'])
            f.write(f'{defect_id},{x:.6f},{y:.6f},{area:.6f},{radius:.6f}\n')


def compute_focus_crop(processed, xc, yc, half_size=FOCUS_HALF_SIZE):
    ny, nx = processed.shape
    x_min = max(0, int(round(xc)) - half_size)
    x_max = min(nx, int(round(xc)) + half_size)
    y_min = max(0, int(round(yc)) - half_size)
    y_max = min(ny, int(round(yc)) + half_size)

    if x_min >= x_max or y_min >= y_max:
        raise RuntimeError(
            f'empty crop for x={xc:.2f}, y={yc:.2f}; check that the coordinates are inside the image'
        )

    focus_image = processed[y_min:y_max, x_min:x_max]
    xc_local = xc - x_min
    yc_local = yc - y_min
    return focus_image, xc_local, yc_local


def save_focus_figure(focus_image, xc_local, yc_local, title, output_basename, vmin=None, vmax=None):
    if vmin is not None and vmax is not None:
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())
    else:
        norm = ImageNormalize(focus_image, stretch=SqrtStretch())

    fig, ax = plt.subplots(figsize=(7.4, 7.6))
    ax.imshow(focus_image, origin='upper', norm=norm)

    circle = Circle((xc_local, yc_local), radius=CIRCLE_RADIUS,
                    edgecolor='red', facecolor='none', linewidth=1.5)
    ax.add_patch(circle)

    ax.set_title(title, fontsize=14, pad=14, wrap=True)
    ax.set_xlim(0, focus_image.shape[1])
    ax.set_ylim(focus_image.shape[0], 0)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if output_basename:
        base = Path(output_basename)
        base.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(base.with_suffix('.png'), dpi=300, bbox_inches='tight')
        fig.savefig(base.with_suffix('.pdf'), bbox_inches='tight')
        print(f' --- saved: {base.with_suffix(".png")}')
        print(f' --- saved: {base.with_suffix(".pdf")}')

    plt.close(fig)


def main():
    args = parse_arguments()

    print(' --- opening input image:', args.input)
    image = tifffile.imread(args.input).astype(float)
    if image.ndim != 2:
        print(f'ERROR: expected a 2D image, got shape {image.shape}', file=sys.stderr)
        return 1

    if args.display_original:
        fig0, ax0 = plt.subplots(figsize=(10, 5))
        ax0.imshow(image, origin='upper')
        ax0.axis('off')
        plt.show()
        plt.close(fig0)

    bkg_estimator = MedianBackground()
    bkg = Background2D(
        image,
        box_size=(50, 50),
        filter_size=(3, 3),
        bkg_estimator=bkg_estimator,
    )
    processed = image - bkg.background

    threshold = 3.0 * bkg.background_rms
    kernel = make_2dgaussian_kernel(3.0, size=5)
    convolved_data = convolve(processed, kernel)

    if args.convolution:
        detection_image = convolved_data
    else:
        detection_image = processed

    finder = SourceFinder(
        npixels=args.npixels,
        deblend=True,
        nlevels=args.nlevels,
        contrast=args.contrast,
        progress_bar=False,
    )

    segment_map = finder(detection_image, threshold)
    if segment_map is None:
        # Still allow pure --focus_area mode.
        tbl = None
        x_col = y_col = None
    else:
        cat = SourceCatalog(processed, segment_map, convolved_data=detection_image)
        tbl = cat.to_table()
        tbl['defect_id'] = np.arange(1, len(tbl) + 1)
        x_col, y_col = get_centroid_column_names(tbl)
        tbl[x_col].info.format = '.2f'
        tbl[y_col].info.format = '.2f'

        segment_area = np.array([to_float(a) for a in cat.segment_area])
        equivalent_radius = np.sqrt(segment_area / np.pi)
        tbl['segment_area_pix2'] = segment_area
        tbl['equivalent_radius_pix'] = equivalent_radius
        tbl['segment_area_pix2'].info.format = '.2f'
        tbl['equivalent_radius_pix'].info.format = '.2f'

        if args.detected_csv:
            save_detected_csv(tbl, x_col, y_col, args.detected_csv)

    if args.display_processed and segment_map is not None:
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
        ax1.imshow(processed, origin='upper', norm=norm)
        ax1.set_title('Background-subtracted Data')
        ax1.axis('off')
        ax2.imshow(segment_map.data, origin='upper', cmap=segment_map.cmap, interpolation='nearest')
        ax2.set_title('Deblended Segmentation Image')
        ax2.axis('off')
        plt.tight_layout()
        plt.show()
        plt.close(fig)

    if args.focus is not None:
        if tbl is None:
            print('ERROR: no sources detected, therefore --focus cannot be used.', file=sys.stderr)
            return 1

        matches = np.where(np.array(tbl['defect_id']) == args.focus)[0]
        if len(matches) == 0:
            print(f'ERROR: selected defect_id {args.focus} is out of range. Detected defect_id values: 1 ... {len(tbl)}', file=sys.stderr)
            return 1
        focus_index = int(matches[0])
        xc = to_float(tbl[x_col][focus_index])
        yc = to_float(tbl[y_col][focus_index])
        print(f' --- focus on defect {args.focus}: x={xc:.2f}, y={yc:.2f}')
        focus_image, xc_local, yc_local = compute_focus_crop(processed, xc, yc)

    elif args.focus_area is not None:
        xc = float(args.focus_area[0])
        yc = float(args.focus_area[1])
        print(f' --- focus on manual area: x={xc:.2f}, y={yc:.2f}, r={CIRCLE_RADIUS:.2f} pix')
        focus_image, xc_local, yc_local = compute_focus_crop(processed, xc, yc)

    else:
        # Detection-only mode.
        return 0

    if args.measure_max:
        print(f'{float(np.nanmax(focus_image)):.12g}')
        return 0

    save_focus_figure(
        focus_image,
        xc_local,
        yc_local,
        args.focus_title,
        args.output_focus,
        vmin=args.vmin,
        vmax=args.vmax,
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())
