#! /usr/bin/env python

"""
Tiff Writing Example - tifffile

This example shows how to use Thorlabs TSI Cameras to write images to a disk using the tifffile library,
see https://pypi.org/project/tifffile/ for more information.

There are many TIFF-writing libraries for python, this example is meant to show how to integrate with tifffile.
The process should generally be the same with most TIFF-writing libraries, but results may vary.

In this example 10 images are going to be taken and saved to a single multipage TIFF file. The program will detect
if the camera has a color filter and will perform color processing if so.

One thing to note is that this program will save TIFFs in the camera's bit depth. Some image viewers may not recognize
this and will show the images as being much darker than expected. If you are experiencing dark images, we recommend
trying out various image viewers designed for scientific imaging such as ThorCam or ImageJ.

"""

try:
    # if on Windows, use the provided setup script to add the DLLs folder to the PATH
    from windows_setup import configure_path
    configure_path()
except ImportError:
    configure_path = None

import os
import tifffile
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="CERNIA")
parser.add_argument("--nframes", type=int, required=True, help="Number of frames to acquire")
parser.add_argument("--exposure", type=float, required=True, help="Exposure in milliseconds")
parser.add_argument("--prefix", type=str, required=True, help="Output file name prefix [.tif]")
parser.add_argument("--average", action="store_true", help="Save average of frames")
parser.add_argument("--no_stack", action="store_true", help="Do not save stack of frames")
parser.add_argument("--pedestal", type=str, required=False, help="Load pedestal image")
args = parser.parse_args()


from thorlabs_tsi_sdk.tl_camera import TLCameraSDK
from thorlabs_tsi_sdk.tl_camera_enums import SENSOR_TYPE

NUMBER_OF_IMAGES = args.nframes
PREFIX = args.prefix
STACK_FILENAME = PREFIX + '.stack.tif'
AVERAGE_FILENAME = PREFIX + '.average.tif'
AVERERR_FILENAME = PREFIX + '.avererr.tif'

TAG_BITDEPTH = 32768
TAG_EXPOSURE = 32769

# delete image if it exists
if os.path.exists(STACK_FILENAME):
    os.remove(STACK_FILENAME)

# load pedestal image
pedestal_image = None
if args.pedestal is not None:
    print(f' --- pedestal image: {args.pedestal}')
    with tifffile.TiffFile(args.pedestal) as tif:
        pedestal_image = tif.asarray()
    
with TLCameraSDK() as sdk:
    cameras = sdk.discover_available_cameras()
    if len(cameras) == 0:
        print("Error: no cameras detected!")

    with sdk.open_camera(cameras[0]) as camera:
        #  setup the camera for continuous acquisition
        camera.frames_per_trigger_zero_for_unlimited = 0
        camera.image_poll_timeout_ms = max(1000, int(args.exposure * 10)) # at least 10x exposure time and not less than 2 s
        camera.arm(2)

        ### set exposure time
        camera.exposure_time_us = int(args.exposure * 1000.)

        ### set LED off
        camera.is_led_on = False
        
        # save these values to place in our custom TIFF tags later
        bit_depth = camera.bit_depth
        exposure = camera.exposure_time_us

        # need to save the image width and height for color processing
        image_width = camera.image_width_pixels
        image_height = camera.image_height_pixels

        # begin acquisition
        camera.issue_software_trigger()
        frames_counted = 0
        H = None
        E = None
        while frames_counted < NUMBER_OF_IMAGES:
            frame = camera.get_pending_frame_or_null()
            if frame is None:
                raise TimeoutError("Timeout was reached while polling for a frame, program will now exit")
            frames_counted += 1
            print(f' --- acquired frame {frames_counted} / {NUMBER_OF_IMAGES}')

            image_data = frame.image_buffer
            image_data = image_data.astype(np.float64)
            if pedestal_image is not None:
                image_data -= pedestal_image
                
            if args.average:
                if H is None:
                    H = image_data
                    E = (image_data * image_data)
                else:
                    H += image_data
                    E += (image_data * image_data)

            if not args.no_stack:
                with tifffile.TiffWriter(STACK_FILENAME, append=True) as tiff:
                    tiff.write(data=image_data.astype(np.float32), # np.int16
                               compression=0,
                               extratags=[(TAG_BITDEPTH, 'I', 1, bit_depth, False),
                                          (TAG_EXPOSURE, 'I', 1, exposure, False)]
                               )

        if args.average:
            W = frames_counted
            h = H / W
            s = np.sqrt(E / W - (h * h))
            e = s / np.sqrt(W)

            h = h.astype(np.float32)
            s = s.astype(np.float32)
            e = e.astype(np.float32)
            
            with tifffile.TiffWriter(AVERAGE_FILENAME) as tiff:
                tiff.write(data=h, dtype=np.float32)
                tiff.write(data=e, dtype=np.float32)
                tiff.write(data=s, dtype=np.float32)
        
        camera.disarm()


