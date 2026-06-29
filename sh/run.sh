#!/bin/bash

echo "Starting"


cd manipulate_images/

python process_image.py \
--input ../files/1originals/run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4_data=diff.tif \
--process remove_column_bias remove_hot_pixels remove_cold_pixels \
--output ../files/2processed/run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4_data=diff_processed.tif

python tif2th2.py \
--input ../files/2processed/run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4_data=diff_processed.tif \
--error ../files/1originals/run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4_data=diffe.tif \
--output ../files/4th2f/run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4_data=diff_processed_th2f.root

cd ../find_centers/
python find_centers.py \
--input ../files/2processed/run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4_data=diff_processed.tif \
--convolution \
--coordinates_root ../files/3coordinates/run=20250826-035846_x==95600_y=50600_z=92600_T=22.0_v=58.4_data=diff_processed_coords.txt

cd ../spot_analysis/
root -l -q 'spot_fit_gaussian_step.C("../files/4th2f/run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4_data=diff_processed_th2f.root","../files/3coordinates/run=20250826-035846_x==95600_y=50600_z=92600_T=22.0_v=58.4_data=diff_processed_coords.txt")'

mv spots.png ../files/fit_results/images_run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4
mv LEGO.png ../files/fit_results/images_run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4
mv fit_results.png ../files/fit_results/images_run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4
mv fit_parameters.png ../files/fit_results/images_run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4
mv fit_results.csv ../files/fit_results/images_run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4

cd ../
echo "Done"


root -l  'spot_luminosity_radius.C("../files/4th2f/run=20250826-035846_x=95600_y=50600_z=92600_T=22.0_v=58.4_data=diff_processed_th2f.root","../files/3coordinates/run=20250826-035846_x==95600_y=50600_z=92600_T=22.0_v=58.4_data=diff_processed_coords.txt")'
