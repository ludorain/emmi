#! /usr/bin/env bash

#while [ $(pgrep -fc ivscan.sh) -gt 0 ]; do
#    echo " there is still at least one running measure_1.sh process "
#    sleep 10
#done

rnohup /emmi/measure/ivscan.sh
