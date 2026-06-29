#! /usr/bin/env python

###
### ~/CERNIA/Scientific_Camera_Interfaces/SDK/Python_Toolkit/thorlabs_tsi_sdk-0.0.8/thorlabs_tsi_sdk/tl_camera.py
###

attributes = '''
 is_led_on
 is_led_supported
 is_cooling_supported
 is_cooling_enabled
 exposure_time_range_us
'''

import os
import tifffile
import argparse
import numpy as np

from thorlabs_tsi_sdk.tl_camera import TLCameraSDK, TLCamera
from thorlabs_tsi_sdk.tl_camera_enums import SENSOR_TYPE

parser = argparse.ArgumentParser(description='CERNIA')
parser.add_argument('--cmd', type=str, required=True, help='Command')
args = parser.parse_args()

camera_class = TLCamera

attribute = None
value = None

cmd_parts = args.cmd.split()
if len(cmd_parts) == 2 and cmd_parts[0] == 'get':
    attribute = cmd_parts[1]
    # check if the attribute exists
    if not hasattr(TLCamera, attribute) or not isinstance(getattr(TLCamera, attribute), property):
        print(f' --- TLCamera: property {attribute} does not exist or is not readable')
        exit(1)
elif len(cmd_parts) == 3 and cmd_parts[0] == 'set':
    attribute = cmd_parts[1]
    try:
        value = int(cmd_parts[2])
    except Exception as e:
        print(e)
        exit(1)
    # check if the attribute has a setter
    prop = getattr(TLCamera, attribute, None)
    if not isinstance(prop, property) or prop.fset is None:
        print(f' --- TLCamera: property {attribute} does not have a setter or does not exist')
        exit(1)
    print(f' --- set {attribute} to {value}')
else:
    print(f' --- invalid command: {cmd}')
    exit(1)

with TLCameraSDK() as sdk:
    cameras = sdk.discover_available_cameras()
    if len(cameras) == 0:
        print("Error: no cameras detected!")

    with sdk.open_camera(cameras[0]) as camera:

        if attribute is not None:
            if value is None:
                try:
                    value = getattr(camera, attribute)
                    print(f' get {attribute}: {value}')
                except Exception as e:
                    print(e)
            else:
                try:
                    setattr(camera, attribute, value)
                    print(f' set {attribute} {value}')
                except Exception as e:
                    print(e)
                    
