Boundary_Condition="PBC"
for Lx in 12 #1.5 2.0 1.0
do
Ly=12
mkdir -p Lx_Ly_${Lx}x${Ly}_${Boundary_Condition}_Run2
cd Lx_Ly_${Lx}x${Ly}_${Boundary_Condition}_Run2
NP=$(echo "4*${Lx}*${Ly}" | bc -l)

mkdir -p NP_${NP}
cd NP_${NP}

for Seed in 1
do
mkdir -p Seed_${Seed}
cd Seed_${Seed}


for U in 4.0 #0.1 0.25 #1.0 3.0 4.0  #18.0 24.0 #1.0 0.5 #10.0 14.0 #0.0 0.5 1.0 1.5 2.0 3.0 4.0 6.0 8.0 10.0 15.0 20.0 #30.0
do
mkdir -p U_${U}
cd U_${U}

for JbyU in 0.25 #0.0 0.1 0.2 0.3
#0.0 0.05 0.1 0.2 0.25 0.3 #0.25
#0.0 0.05 0.1 0.2 0.25 0.3
do
mkdir -p JbyU_${JbyU}
cd JbyU_${JbyU}

J=$(echo "${U}*${JbyU}" | bc -l)
J=$(printf "%1.5f" ${J})

UPrime=$(echo "${U} - 2.0*${J}" | bc -l)
UPrime=$(printf "%1.5f" ${UPrime})

#16 2.1 2.11 2.12 2.13 2.14 2.15 2.18 2.2 2.22 2.25 2.28 2.3 2.32 2.35 2.38
#8 2.1 2.12 2.14 2.15 2.18 2.2 2.22 2.25 2.28 2.26 2.27 2.28 2.29 2.3 2.32 2.35 2.38 2.4
#24 2.1 2.12 2.14 2.15 2.18 2.2 2.22 2.25 2.28 2.3 2.32 2.35 2.38
for lambda in 3.0 #2.22 2.25 2.28 2.32 2.355 2.36 2.365 2.375 2.11 2.12 2.13 2.14 2.15 2.18 2.2 2.3 2.35 2.37 2.38
      	#2.355 2.36 2.365 2.37 2.375 #2.1 2.11 2.13 2.14 2.16 #2.12 2.15 2.18 2.2 2.22 2.25 2.28 2.3 2.32 2.35 2.38
	#2.1 2.11 2.12 2.13 2.14 2.15 2.18 2.2 2.22 2.25 2.28 2.3 2.32 2.35 2.38 #2.11 2.13 2.14 #2.26 2.27 2.29 #2.12 2.18 2.22 2.28 2.32 2.38   #2.1 2.12 2.15  2.18 2.2 2.22 2.25 2.28 2.3 2.32 2.35 2.38   #0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2.0 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5 2.55 2.6 2.8 #0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 2.0 

#For U=0.5, 1.0
#0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2.0 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5 2.55 2.6 2.8 

do
mkdir -p lambda_${lambda}
cd lambda_${lambda}

for dis_str in 0.25 0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.5 5.75 6.25 6.5 6.75 7.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 6.0 
do
mkdir -p dis_str_${dis_str}
cd dis_str_${dis_str}


for dis_seed in 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 #21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40  #1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
mkdir -p dis_seed_${dis_seed}
cd dis_seed_${dis_seed}



input=input_run.inp
rm ${input}
cp ../../../../../../../../input_template_${Boundary_Condition}_Run2.inp ${input}
cp ../../../../../../../../SelfConsistent .

#chmod 777 SelfConsistent
#cp ../../../../../Runs/J_${J_exchange}/Dis_${Dis}/unbiased_OP_seed_${unbiased_OP_seed}/Disorder_seed_${Disorder_seed}/output_Local_* .


########################################
sed -i -e "s/VALUE_LX/${Lx}/g" $input
sed -i -e "s/VALUE_LY/${Ly}/g" $input
sed -i -e "s/VALUE_NP/${NP}/g" $input
sed -i -e "s/VALUE_SEED/${Seed}/g" $input
sed -i -e "s/VALUE_HUND/${J}/g" $input
sed -i -e "s/VALUE_U/${U}/g" $input
sed -i -e "s/VALUE_PRIME_U/${UPrime}/g" $input
sed -i -e "s/VALUE_LAMBDA/${lambda}/g" $input
sed -i -e "s/VALUESTRENGTH_DISORDER/${dis_str}/g" $input
sed -i -e "s/VALUESEED_DISORDER/${dis_seed}/g" $input
########################################

sub="job_GS_${Seed}_U${U}"
rm $sub
#     $sub
        echo "#!/bin/sh"                        >> $sub
        echo "#PBS -N RUN_HF"                            >>$sub
        echo "#PBS -l nodes=1:ppn=4"                   >> $sub
        echo "#PBS -l walltime=100:50:00"               >> $sub
        echo "#PBS -l vmem=44gb" >> $sub
        echo "#PBS -q batch"  >> $sub
        echo "hostname"         >> $sub
        echo "#PBS -r n"                            >>$sub
        echo 'cd $PBS_O_WORKDIR'                        >> $sub
        echo "date"                             >> $sub
        echo "time ./SelfConsistent ${input} > out_run.txt"              >> $sub
        echo "date"                             >> $sub

##################### COMPILE THE CODE #########################################
qsub $sub

#echo "$sub submitted "


echo "Lx=${Lx} Ly=${Ly} NP=${NP} Seed_OP=${Seed} U=${U} Lambda=${lambda} Dis_STR=${dis_str} Dis_Seed=${dis_seed} submitted"
#printf "%1.5f" ${J}
#echo ""


cd ..
done #Dis_seed

cd ..
done #Dis_Str

cd ..
done #lambda

cd ..  #JbyU
done

cd ..  #U
done

cd ..  #Seed
done

cd ..  #NP

cd .. #L
done
