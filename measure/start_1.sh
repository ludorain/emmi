#! /usr/bin/env bash

export EMMI_BOARD="33"      # Number of the SiPM carrier board
export EMMI_SENSOR="B4"     # Position of the sensor in the SiPM matrix
export EMMI_TEMPS="20"      # Temperature
export EMMI_VOVERS="7"      # Overvoltage
export EMMI_POSITIONS=ALL   # ALL (all sensor) or GOLD (only one image)
export EMMI_LIGHTONLY=false # false=performs all measurements, true=optical image only
export EMMI_RMRAW=true      # true=removes all raw images, false=saves all images 

while [ $(pgrep -fc measure_1.sh) -gt 0 ]; do
    echo " there is still at least one running measure_1.sh process "
    sleep 10
done

rnohup /emmi/measure/measure_1.sh
