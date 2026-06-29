import skimage
import matplotlib.pyplot as plt

def read_image(filename, verbose=False):
    print(' --- read image:', filename)
    image = skimage.io.imread(filename, as_gray=True)
    if verbose: 
        print('     object type: ', type(image))
        print('     image shape:', image.shape)
        print('     image dtype:', image.dtype)
    return image

def save_image(filename, image):
    print(' --- save image:', filename)
    image = skimage.io.imsave(filename, image)

def print_png(filename, image, cmap, vmin, vmax):
    print(' --- print png:', filename)
    plt.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.axis('off')
    plt.savefig(filename, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()

