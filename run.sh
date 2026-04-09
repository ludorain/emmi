#!/bin/bash

echo "Starting"

cd find_centers/
python3 find_centers.py --input 
cd ../manipulate_images/
python3 process-image.py
cd ../spot_analysis/
root -l -q 'spot_fit_gaussian_circular_noB.C("results/fit_results.csv")'
cd ../
echo "Done"