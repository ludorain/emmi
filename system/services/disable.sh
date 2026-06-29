#! /usr/bin/env bash

systemctl disable emmi.target

### LM73 services
systemctl disable emmi-lm73-read.service
systemctl disable emmi-lm73-influxdb.service
