#! /usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np

import sys
import os

current_dir = os.path.dirname(__file__)
libs_path = os.path.abspath(os.path.join(current_dir, '..', 'libs'))
sys.path.insert(0, libs_path)
import emmi

def display(image, angle, dist):

    plt.imshow(image, cmap='gray')
    plt.axis('off')

    (x0, y0) = dist * np.array([np.cos(angle), np.sin(angle)])
    plt.axline((x0, y0), slope=np.tan(angle + np.pi / 2), color='red', linestyle=(0, (5, 5)))

    plt.show()    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="D")
    parser.add_argument("--input", required=True, help="Input TIF file name")
    parser.add_argument("--display", action="store_true", help="Display process")    
    args = parser.parse_args()

    ### read image
    image = emmi.io.read_image(args.input)

    ### get rotation angle
    rotation = emmi.rotation.measure_rotation(image)
    print(f' --- rotation angle: {rotation["correction_degrees"]} deg')

    ### display process
    if args.display:
        display(image, rotation['angle'], rotation['distance'])