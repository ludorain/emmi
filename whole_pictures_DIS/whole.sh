cd ..
cd find_centers/
python find_centers.py \
--input ../whole_pictures_DIS/originals/20250823-232851.stitched.data=diff-denoised.tif \
--convolution \
--display_original \
--circles \
--coordinates_root ../whole_pictures_DIS/result_count/0250823-232851.stitched.data=diff-denoised.txt

python find_centers.py \
--input ../whole_pictures_DIS/originals/20250927-031407.stitched.data=diff-denoised.tif \
--convolution \
--display_original \
--circles \
--coordinates_root ../whole_pictures_DIS/result_count/20250927-031407.stitched.data=diff-denoised.txt

python find_centers.py \
--input ../whole_pictures_DIS/originals/20260512-123812.stitched.data=diff-denoised.tif \
--convolution \
--display_original \
--circles \
--coordinates_root ../whole_pictures_DIS/result_count/20260512-123812.stitched.data=diff-denoised.txt
cd ..