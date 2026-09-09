import numpy as np
import skimage


def measure_rotation(image):
    thresh = skimage.filters.threshold_otsu(image)
    image = image > thresh
    image = skimage.filters.sobel(image)
    tested_angles = np.linspace(-np.pi / 100, np.pi / 100, 500, endpoint=False)
    h, theta, d = skimage.transform.hough_line(image, theta=tested_angles)
    lines = list(zip(*skimage.transform.hough_line_peaks(h, theta, d)))
    angle, distance = lines[0][1:]
    degrees = np.rad2deg(angle)
    correction_degrees = degrees - 90 if degrees > 45 else degrees
    return {
        'angle': angle,
        'distance': distance,
        'correction_degrees': correction_degrees,
    }


def get_rotation_angle(image):
    return measure_rotation(image)['correction_degrees']


def rotate_image(image, degrees):
    print(f' --- rotate image: {degrees} degrees')
    rotated = skimage.transform.rotate(image, degrees, resize=False)
    return rotated