import skimage
import numpy as np

save_intensity_range = None

def fix_underflows(image, thresh = -1.e3):
    image = np.where(image < thresh, np.max(image), image)
    return image

def contrast_stretching(image):
    print(' --- contrast_stretching')
    pmin, pmax = np.percentile(image, (1, 99))
    image = skimage.exposure.rescale_intensity(image, in_range=(pmin, pmax))
    return image

def rescale_intensity(image):
    print(' --- rescale_intensity')
    global save_intensity_range
    save_intensity_range = (np.min(image), np.max(image))
    image = skimage.exposure.rescale_intensity(image, in_range='image', out_range=(0, 1))
    return image

def restore_intensity(image):
    print(' --- restore_intensity')
    global save_intensity_range
    if save_intensity_range is not None:
        image = skimage.exposure.rescale_intensity(image, in_range='image', out_range=save_intensity_range)
        save_intensity_range = None
    return image    
    
def equalize_adapthist(image):
    print(' --- equalize_adapthist')
    save_intensity_range = (np.min(image), np.max(image))
    image = skimage.exposure.rescale_intensity(image, in_range='image', out_range=(0, 1))
    image = skimage.exposure.equalize_adapthist(image)
    image = skimage.exposure.rescale_intensity(image, in_range='image', out_range=save_intensity_range)
    return image

def equalize_hist(image):
    print(' --- equalize_hist')
    image = skimage.exposure.equalize_hist(image)
    return image

def threshold_local(image):
    print(' --- threshold_local')
    thresh = skimage.filters.threshold_local(image, block_size=35)
    image = image > thresh
    return image

def threshold_otsu(image):
    print(' --- threshold_otsu')
    thresh = skimage.filters.threshold_otsu(image)
    image = image > thresh
    return image

def sobel(image):
    print(' --- sobel')
    image = skimage.filters.sobel(image)
    return image

def canny(image):
    print(' --- canny')
    image = skimage.feature.canny(image)
    return image

def remove_column_bias(image):
    print(' --- remove columns bias')
    column_bias = np.median(image, axis=0);
    return (image - column_bias)

def remove_hot_pixels(image):
    print(' --- remove hot pixels')
    local_median = skimage.filters.median(image, skimage.morphology.disk(1))
    # identify and replace hot pixels
    sigma = skimage.restoration.estimate_sigma(image)
    hot_pixels = (image - local_median) > (5. * sigma)
    image_out = image.copy()
    image_out[hot_pixels] = local_median[hot_pixels]
    return image_out

def remove_cold_pixels(image):
    print(' --- remove cold pixels')
    local_median = skimage.filters.median(image, skimage.morphology.disk(1))
    # identify and replace hot pixels
    sigma = skimage.restoration.estimate_sigma(image)
    cold_pixels = (image - local_median) < (5. * sigma)
    image_out = image.copy()
    image_out[cold_pixels] = local_median[cold_pixels]
    return image_out

def denoise_nl_means(image):
    print(' --- denoise nl means')
    sigma = skimage.restoration.estimate_sigma(image)
    # apply Non-Local Means Denoising
    return skimage.restoration.denoise_nl_means(image,
                                                h=1.15 * sigma,
                                                fast_mode=True,
                                                patch_size=5,
                                                patch_distance=6,
                                                channel_axis=None)
