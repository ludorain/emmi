#!/bin/bash

tools_light_on/measure-rotation.py --input /Users/ludovicarainero/emmi/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/1originals/run=20260513-032137_x=0_y=0_z=0_data=light.tif --display

tools_light_on/rotate-image.py --input /Users/ludovicarainero/emmi/DATA_irradiated/before_annealing/A1_T=20_run=20260513-032137/1originals/run=20260513-032137_x=0_y=0_z=0_data=light.tif --angle 0.3384 --output /Users/ludovicarainero/emmi/prova_alignment/A1_T=20_run=20260513-032137-straight.tif

#immagine post annealing
tools_light_on/measure-rotation.py --input /Users/ludovicarainero/emmi/DATA_irradiated/annealing_T=75_h=5/A1_T=20_run=20260520-085312/1originals/run=20260520-085312_x=0_y=0_z=0_data=light.tif --display

tools_light_on/rotate-image.py --input /Users/ludovicarainero/emmi/DATA_irradiated/annealing_T=75_h=5/A1_T=20_run=20260520-085312/1originals/run=20260520-085312_x=0_y=0_z=0_data=light.tif --angle 0.468 --output /Users/ludovicarainero/emmi/prova_alignment/A1_T=20_run=20260520-085312-straight-moving.tif

#misura la distanza tra le due immagini

tools_light_on/measure-shift.py --input /Users/ludovicarainero/emmi/prova_alignment/A1_T=20_run=20260513-032137-straight-reference.tif /Users/ludovicarainero/emmi/prova_alignment/A1_T=20_run=20260520-085312-straight-moving.tif

tools_light_on/shift-image.py --input /Users/ludovicarainero/emmi/prova_alignment/A1_T=20_run=20260520-085312-straight-moving.tif --shift 1. -4. --output /Users/ludovicarainero/emmi/prova_alignment/A1_T=20_run=20260520-085312-straight-moving-shifted.tif