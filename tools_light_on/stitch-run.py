#! /usr/bin/env python

import emmi
import argparse
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True, help='Wildcard pattern (e.g. \'./data/*.txt\')')
    parser.add_argument('--output', type=str, required=False, help='Output file tagname')
    parser.add_argument('--tags', type=str, nargs='+', required=True, help='Data tags to process')
    parser.add_argument('--display', action='store_true', help='Display stitching process')
    parser.add_argument('--final', action='store_true', help='Display final image')
    return parser.parse_args()

def stitch_display(image1, image2, stitched, cmap=None, vmin=None, vmax=None):

    fig = plt.figure(figsize=(10, 5), layout="constrained")
    gs = GridSpec(2, 2, figure=fig)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.imshow(emmi.process.contrast_stretching(image1), cmap=cmap, vmin=vmin, vmax=vmax)
    ax1.set_title('image 1')
    ax1.axis('off')
    
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.imshow(emmi.process.contrast_stretching(image2), cmap=cmap, vmin=vmin, vmax=vmax)
    ax2.set_title('image 2')
    ax2.axis('off')

    if stitched is not None:
        ax3 = fig.add_subplot(gs[1, :])
        ax3.imshow(emmi.process.contrast_stretching(stitched), cmap=cmap, vmin=vmin, vmax=vmax)
        ax3.set_title('stitched')
        ax3.axis('off')
    
    plt.tight_layout()
    plt.show()    

def crop_images(image1, image2, seed):
    if seed[1] > 0:
        image1 = image1[ : , int(seed[1]):]
        image2 = image2[ : , :-int(seed[1])]
    elif seed[1] < 0:
        image2 = image2[ : , -int(seed[1]):]
        image1 = image1[ : , :int(seed[1])]
        
    if seed[0] > 0:
        image1 = image1[int(seed[0]): , :]
        image2 = image2[:-int(seed[0]) , :]
    elif seed[0] < 0:
        image2 = image2[-int(seed[0]): , :]
        image1 = image1[:int(seed[0]) , :]

    return image1, image2

def measure_shift(first, second):
    pos1, image1 = first
    pos2, image2 = second
#    eshift = (pos1 - pos2) * 844. / 5000.
    exshift = int ( (pos1[0] - pos2[0]) * 844. / 5000. )
    eyshift = int ( (pos1[1] - pos2[1]) * 844. / 5000. )
    eshift = (exshift, eyshift)
    image1, image2 = crop_images(image1, image2, eshift)
    shift = eshift + emmi.stitching.phase_correlate(image1, image2)
    return shift

def measure_shifts(image_data):
    shifts = {}
    if len(image_data) < 2:
        print(' --- measure_shifts requires at least two images')
        return None
    for first, second in zip(image_data, image_data[1:]):
        pos1, image1 = first
        pos2, image2 = second
        shift = measure_shift(first, second)
        shifts[(pos1, pos2)] = shift
    return shifts

def stitch_images(image_data, shifts, display, panorama):
    if len(image_data) < 2:
        print(' --- stitch_images requires at least two images')
        return None
    for first, second in zip(image_data, image_data[1:]):
        pos1, image1 = first
        pos2, image2 = second
        shift = shifts[(pos1, pos2)]
        stitched = emmi.stitching.stitch_images(image1, image2, shift, panorama=panorama)
        first[1] = stitched
        if display:
            stitch_display(image1, image2, stitched)
    image_data = image_data[:-1]
    return image_data

if __name__ == "__main__":
    args = parse_arguments()

    # build database and coordinates
    database = emmi.database.build_database(args.input)
    coords = emmi.database.build_coordinates(database)

    tags = {'light'}
    tags.update(args.tags)
    
    ### first create stitched rows

    stitched_row_data = {}
    for tag in tags:
        stitched_row_data[tag] = []

    # loop over rows
    for row in coords:
        
        # load row images
        row_data = {}
        for tag in tags:
            row_data[tag] = []
            for pos in row:
                y, x = pos
                conditions = { 'x':x , 'y':y , 'data':tag }
                filename = emmi.database.get_filename(database, conditions)
                image = emmi.io.read_image(filename)
                row_data[tag].append( [pos, image] )

        # loop until we are left with only one image in the list
        while len(row_data['light']) > 1:
            # measure shifts with light data
            shifts = measure_shifts(row_data['light'])
            # stitch images of all tags
            for tag in tags:
                row_data[tag] = stitch_images(row_data[tag], shifts, display=args.display, panorama=(True,False))

        # populate list of stitched row data
        for tag in tags:
            stitched_row_data[tag].append(row_data[tag][0])

    ### then stitch the rows togheter

    # loop until we are left with only one image in the list
    while len(stitched_row_data['light']) > 1:
        # measure shifts with light
        shifts = measure_shifts(stitched_row_data['light'])        
        # stitch images of all tags
        for tag in tags:
            stitched_row_data[tag] = stitch_images(stitched_row_data[tag], shifts, display=args.display, panorama=(False,True))

    # retrieve stitched data and images
    stitched_data = {}
    stitched_image = {}
    for tag in tags:
        stitched_data[tag] = stitched_row_data[tag][0]
        stitched_image[tag] = stitched_data[tag][1]
        if args.output is not None:
            outfilename = args.output + '.data=' + tag + '.tif'
            emmi.io.save_image(outfilename, stitched_image[tag])
        if args.final:
            fig = plt.figure(figsize=(10, 10), layout="constrained")
            plt.imshow(emmi.process.contrast_stretching(stitched_image[tag]))
            plt.title(tag)
            plt.axis('off')
            plt.tight_layout()
            plt.show()

