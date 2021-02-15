

theta=4.0
ref_site=9

g++ Gradientized_spins_triangular.cpp 

for eps in 1.00 2.00 2.50 3.00 3.50 4.00 4.50 5.00 5.50 6.00 6.50 7.00 7.50 8.00 8.50 9.00 9.50 10.00 10.50 11.00 12.00 13.00 14.00 15.00 18.00 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 70.0 80.0 90.0 100.0
do

file=SS_refsite${ref_site}_theta${theta}_eps${eps}

./a.out 4 4 ${file} ${eps}
pdflatex ${file}.tex

done
