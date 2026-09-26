cd find_centers/
python find_centers.py \
--input ../whole_pictures/originals/20250823-232851.stitched.data=diff-denoised.tif \
--convolution \
--display_original \
--circles \
--coordinates_root ../whole_pictures/result_count/0250823-232851.stitched.data=diff-denoised.txt

python find_centers.py \
--input ../whole_pictures/originals/20250927-031407.stitched.data=diff-denoised.tif \
--convolution \
--display_original \
--circles \
--coordinates_root ../whole_pictures/result_count/20250927-031407.stitched.data=diff-denoised.txt

python find_centers.py \
--input ../whole_pictures/originals/stitched.data=diff-denoised.tif \
--convolution \
--display_original \
--circles \
--coordinates_root ../whole_pictures/result_count/20250927-031407.stitched.data=diff-denoised.txt
cd ..