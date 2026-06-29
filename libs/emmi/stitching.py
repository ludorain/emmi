import numpy as np
import skimage
import scipy
from emmi import process

def phase_correlate(image1, image2, reference_mask=None, moving_mask=None):
    print(' --- phase correlate images ')
    h = min(image1.shape[0], image2.shape[0])
    w = min(image1.shape[1], image2.shape[1])
    image1 = image1[:h, :w]
    image2 = image2[:h, :w]

    image1 = process.fix_underflows(image1)
    image2 = process.fix_underflows(image2)
    image1 = process.equalize_adapthist(image1)
    image2 = process.equalize_adapthist(image2)
    if reference_mask is not None:
        reference_mask = reference_mask[:h, :w]
    if moving_mask is not None:
        moving_mask = moving_mask[:h, :w]
#    image1 = skimage.feature.canny(image1, sigma=1)
#    image2 = skimage.feature.canny(image2, sigma=1)
    image1 = skimage.filters.sobel(image1)
    image2 = skimage.filters.sobel(image2)
    shift, error, diffphase = skimage.registration.phase_cross_correlation(image1, image2, upsample_factor=1,
                                                                           reference_mask=reference_mask, moving_mask=moving_mask)
    print(f'     detected subpixel offset (y, x): {shift}')
    return shift

def shift_image(image, shift):
    print(f' --- shift image with offset (y, x): {shift}')
    image = scipy.ndimage.shift(image, shift=shift, order=1, mode='constant', cval=0.0)
    return image

def stitch_images(image1, image2, shift=None, panorama=(False,False)):
    print(f' --- stitch images with offset (y, x): {shift}')
    if shift is None:
        shift = phase_correlate(image1, image2)

    dy, dx = shift
    h1, w1 = image1.shape
    h2, w2 = image2.shape
    
    # Separate integer and fractional parts
    dx_int, dy_int = int(np.floor(dx)), int(np.floor(dy))
    dx_frac, dy_frac = dx - dx_int, dy - dy_int
    
    # Fractional shift with trimming
    image2 = scipy.ndimage.shift(image2, shift=(dy_frac, dx_frac), order=1, mode='constant', cval=0.0)
#    image2 = image2[1:-1, 1:-1]   # trim unreliable borders
#    h2, w2 = image2.shape
    
    # Bounding box
    x_min = min(0, dx_int)
    y_min = min(0, dy_int)
    x_max = max(w1, dx_int + w2)
    y_max = max(h1, dy_int + h2)
    
    canvas_w, canvas_h = x_max - x_min, y_max - y_min
    canvas = np.zeros((canvas_h, canvas_w), dtype=np.float32)
    weight_map = np.zeros((canvas_h, canvas_w), dtype=np.float32)
    
    # Place image1
    x1, y1 = -x_min, -y_min
    canvas[y1:y1+h1, x1:x1+w1] += image1
    weight_map[y1:y1+h1, x1:x1+w1] += 1.0
    
    # Place shifted image2
    x2, y2 = x1 + dx_int, y1 + dy_int
    canvas[y2:y2+h2, x2:x2+w2] += image2
    weight_map[y2:y2+h2, x2:x2+w2] += 1.0
    
    # Average overlaps
    mask = weight_map > 0
    canvas[mask] /= weight_map[mask]
    
    if panorama[0]:
        canvas = canvas[ abs(dy_int):canvas.shape[0]-abs(dy_int), : ]
    if panorama[1]:
        canvas = canvas[ : , abs(dx_int):canvas.shape[1]-abs(dx_int)]
    
    print('     stitched image shape:', canvas.shape)
    return canvas
    

def stitch_images_old(image1, image2, shift=None, panorama=(False,False)):
    print(f' --- stitch images with offset (y, x): {shift}')
    if shift is None:
        shift = phase_correlate(image1, image2)
    dy, dx = shift
    h1, w1 = image1.shape
    h2, w2 = image2.shape
    # Separate integer and fractional parts for img2 shift
    dx_int = int(np.floor(dx))
    dy_int = int(np.floor(dy))
    dx_frac = dx - dx_int
    dy_frac = dy - dy_int
    # Apply fractional shift to img2 (subpixel alignment)
    image2 = scipy.ndimage.shift(image2, shift=(dy_frac, dx_frac), order=1, mode='constant', cval=0.0)
    # Determine bounding box
    x_min = min(0, dx_int)
    y_min = min(0, dy_int)
    x_max = max(w1, dx_int + w2)
    y_max = max(h1, dy_int + h2)
    canvas_w = int(x_max - x_min)
    canvas_h = int(y_max - y_min)
    # Create float32 canvas and weight map
    canvas = np.zeros((canvas_h, canvas_w), dtype=np.float32)
    weight_map = np.zeros((canvas_h, canvas_w), dtype=np.float32)
    # Place img1
    x1 = -x_min
    y1 = -y_min
    canvas[y1:y1+h1, x1:x1+w1] += image1
    weight_map[y1:y1+h1, x1:x1+w1] += 1.
    # Place shifted img2
    x2 = x1 + dx_int
    y2 = y1 + dy_int
    canvas[y2:y2+h2, x2:x2+w2] += image2
    weight_map[y2:y2+h2, x2:x2+w2] += 1.0
    # Average overlaps
    mask = weight_map > 0
    canvas[mask] /= weight_map[mask]

    if panorama[0]:
        canvas = canvas[ abs(dy_int):canvas.shape[0]-abs(dy_int), : ]
    if panorama[1]:
        canvas = canvas[ : , abs(dx_int):canvas.shape[1]-abs(dx_int)]
    
    print('     stitched image shape:', canvas.shape)
    return canvas
