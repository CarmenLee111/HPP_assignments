
make clean
make


# echo 4 body result
# ./galsim 4 input_data/circles_N_4.gal 1 0.00001 0 0
# ./galsim 2 input_data/circles_N_2.gal 1 0.00001 0 0

# export NUM_THREADS=4

# ./galsim 4 input_data/circles_N_4.gal 1 0.00001 0 0 1
# ./galsim 4 input_data/circles_N_4.gal 1 0.00001 0 0 1


# Checking accuracy

echo Compare 10 bodies results
./galsim 10 input_data/ellipse_N_00010.gal 200 0.00001 0 0 1
./compare_gal_files/comp 10 ref_output_data/ellipse_N_00010_after200steps.gal result.gal
rm result.gal  

for i in 1 2 3 4
do 
    # export OMP_NUM_THREADS=$i
    echo Timing 5000 bodies 100 steps with $i threads
    time ./galsim 5000 input_data/ellipse_N_05000.gal 100 0.00001 0.25 0 $i
done
 
# echo Compare 100 bodies results
# ./galsim 100 input_data/ellipse_N_00100.gal 200 0.00001 0 0 0
# ./compare_gal_files/comp 100 ref_output_data/ellipse_N_00100_after200steps.gal result.gal
# rm result.gal



# echo Compare 100 bodies results
# ./galsim 100 input_data/ellipse_N_00100.gal 200 0.00001 0 0 4
# ./compare_gal_files/comp 100 ref_output_data/ellipse_N_00100_after200steps.gal result.gal
# rm result.gal

# echo Compare 1000 bodies results
# ./galsim 1000 input_data/ellipse_N_01000.gal 200 0.00001 0.25 0 4
# ./compare_gal_files/comp 1000 ref_output_data/ellipse_N_01000_after200steps.gal result.gal
# rm result.gal

# ./galsim 02000 input_data/ellipse_N_02000.gal 800 0.00001 0.25 1 4

# echo test Result

# touch testResult.txt
# echo Barne-Hut Complexity with parallelization | tee -a testResult.txt


# for i in 1 2 4 8 
# do 
#     echo 2000 bodies 200 steps with $i threads #| tee -a testResult.txt
#     (time ./galsim 2000 input_data/ellipse_N_02000.gal 200 0.00001 0.25 0 $i) 2>> testResult.txt
#     rm result.gal
# done


# for i in 1 2 4 8
# do
#     echo Timing 3000 bodies 100 steps with $i threads # | tee -a testResult.txt
#     (time ./galsim 3000 input_data/ellipse_N_03000.gal 100 0.00001 0.25 0 $i) 2>> testResult.txt
#     echo ---------------------- >> testResult.txt
# done

# for i in 1 2 4 8 16 32
# do
#     echo  
#     echo Timing 5000 bodies 100 steps with $i threads
#     time ./galsim 5000 input_data/ellipse_N_05000.gal 100 0.00001 0.25 0 $i
# done



# for i in 0 1 2 
# do 
#     echo 1500 bodies results | tee -a testResult.txt
#     (time ./galsim 1500 input_data/ellipse_N_01500.gal 200 0.00001 0.25 0) 2>> testResult.txt
#     rm result.gal
# done

# for i in 0 1 2 
# do 
#     echo 2000 bodies results | tee -a testResult.txt
#     (time ./galsim 2000 input_data/ellipse_N_02000.gal 200 0.00001 0.25 0 8) 2>> testResult.txt
#     rm result.gal
# done

# for i in 0 1 2 
# do 
#     echo 3000 bodies results | tee -a testResult.txt
#     (time ./galsim 3000 input_data/ellipse_N_03000.gal 200 0.00001 0.25 0) 2>> testResult.txt
#     rm result.gal
# done






