

theta=4.0
ref_site=0

g++ Gradientized_spins_triangular.cpp 

for eps in 2.0 4.5 6.5 7.0 7.5 11.0 13.0 20.0 #2.5 3.5 4.5 5.0 5.5 6.5 7.0 7.5 8.5 9.0 9.5 10.5 11.0 12.0 13.0 14.0 1.0 2.0 4.0 6.0 8.0 10.0 15.0 18.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 70.0 80.0 90.0 100.0
do

file=denden_refsite${ref_site}_theta${theta}_eps${eps}

./a.out 4 4 ${file} ${eps}
pdflatex ${file}.tex

done
