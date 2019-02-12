
make clean
make

echo Compare 10 bodies results

./galsim 10 input_data/ellipse_N_00010.gal 200 0.00001 0
./compare_gal_files/comp 10 ref_output_data/ellipse_N_00010_after200steps.gal result.gal
rm result.gal

echo Compare 100 bodies results
./galsim 100 input_data/ellipse_N_00100.gal 200 0.00001 0
./compare_gal_files/comp 100 ref_output_data/ellipse_N_00100_after200steps.gal result.gal
rm result.gal

echo Compare 100 bodies results
./galsim 100 input_data/ellipse_N_00100.gal 200 0.00001 0
./compare_gal_files/comp 100 ref_output_data/ellipse_N_00100_after200steps.gal result.gal
rm result.gal

echo Compare 500 bodies results
./galsim 500 input_data/ellipse_N_00500.gal 200 0.00001 0
./compare_gal_files/comp 500 ref_output_data/ellipse_N_00500_after200steps.gal result.gal
rm result.gal

echo Compare 1000 bodies results
./galsim 1000 input_data/ellipse_N_01000.gal 200 0.00001 0
./compare_gal_files/comp 1000 ref_output_data/ellipse_N_01000_after200steps.gal result.gal
rm result.gal

echo Compare 2000 bodies results
./galsim 2000 input_data/ellipse_N_02000.gal 200 0.00001 0
./compare_gal_files/comp 2000 ref_output_data/ellipse_N_02000_after200steps.gal result.gal
rm result.gal

echo Compare 3000 bodies results
./galsim 3000 input_data/ellipse_N_03000.gal 100 0.00001 0
./compare_gal_files/comp 3000 ref_output_data/ellipse_N_03000_after100steps.gal result.gal
rm result.gal


touch testResult.txt
for i in 0 1 2 3
do
    echo Timing no optimization 10 bodies 200 steps >> testResult.txt
    (time ./galsim 10 input_data/ellipse_N_00010.gal 200 0.00001 0) 2>> testResult.txt
    echo ---------------------- >> testResult.txt 
done


