cd ..
cd ..

cd find_centers/
python find_defects_NOisolated_fixedR.py \
--input /Users/ludovicarainero/emmi/new_device/whole_device/originals/A1_total_20250823-232851.stitched.data=diff-denoised.tif \
--device "A1 new device" \
--circles \
--save \
--save_dir ../new_device/whole_device/A1_images \
--coordinates_root ../new_device/whole_device/result_count/A1_total_20250823-232851.stitched.data=diff-denoised.txt
