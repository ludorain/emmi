import numpy as np

from . import io
from . import database
from . import process
from . import stitching

def get_rotation_angle(image):
    # detect edges
#    image = local_thresholding(image)
    thresh = skimage.filters.threshold_otsu(image)
    image = image > thresh
    image = skimage.filters.sobel(image)
    # straight-line Hough transform
    tested_angles = np.linspace(-np.pi / 100, np.pi / 100, 500, endpoint=False)
    h, theta, d = skimage.transform.hough_line(image, theta=tested_angles)
    lines = list(zip(*skimage.transform.hough_line_peaks(h, theta, d)))
    return lines[0][1:]

def rotate_image(image, degrees):
    print(f' --- rotate image: {degrees} degrees')
    rotated = skimage.transform.rotate(image, degrees, resize=False)
    return rotated

