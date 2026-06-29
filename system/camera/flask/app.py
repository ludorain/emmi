#! /usr/bin/env python

from thorlabs_tsi_sdk.tl_camera import TLCameraSDK
from thorlabs_tsi_sdk.tl_camera_enums import SENSOR_TYPE

camera = None
exposure_time_ms = 300

from flask import Flask, jsonify, render_template
from flask_cors import CORS
from flask_socketio import SocketIO

app = Flask(__name__)
CORS(app)  # Enable CORS for frontend communication
socketio = SocketIO(app, cors_allowed_origins="*")

import base64
import numpy as ns
import cv2

@app.route('/')
def index():
    return render_template('index.html')

@socketio.on("arm")
def arm():
    print(' --- arm ')
    camera.arm(2)
 
@socketio.on("disarm")
def disarm():
    print(' --- disarm ')
    camera.disarm()

@socketio.on("trigger")
def disarm():
    print(' --- trigger ')
    camera.issue_software_trigger()

running = False
    
@socketio.on("acquire")
def acquire():
    frame = camera.get_pending_frame_or_null()
    if frame is None:
        print(' --- timeout ')
        return
    print(' --- acquired image')
    image = frame.image_buffer
    _, buffer = cv2.imencode('.png', image)
    img_base64 = base64.b64encode(buffer).decode('utf-8')
    socketio.emit('image_data', { 'image': img_base64 })
   
if __name__ == "__main__":

    with TLCameraSDK() as sdk:
        cameras = sdk.discover_available_cameras()
        if len(cameras) == 0:
            print("Error: no cameras detected!")
            exit(1)

        with sdk.open_camera(cameras[0]) as c:
            camera = c
            
            #  setup the camera for continuous acquisition
            camera.frames_per_trigger_zero_for_unlimited = 0
            camera.image_poll_timeout_ms = max(1000, int(exposure_time_ms * 10)) # at least 10x exposure time and not less than 2 s
            
            ### set exposure time
            camera.exposure_time_us = int(exposure_time_ms * 1000)
            
            ### set LED off
            camera.is_led_on = False

            socketio.run(app, debug=True, host="0.0.0.0", port=5001, use_reloader=False)
            
