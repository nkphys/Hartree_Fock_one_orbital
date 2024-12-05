#!/bin/bash
#SBATCH --job-name=Run_w_l
#SBATCH -A st-ianaffle-1
#SBATCH -t 23:00:00
#SBATCH --output=output.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
cd $SLURM_SUBMIT_DIR
module load gcc/9.4.0
module load intel-mkl/2020.4.304
date
time ./k_space_SelfConsistency Generic2dLattice input_run.inp > out_run.txt
rm Nw0.0001000000.txt
rm k_space_SelfConsistency
rm Bands*
rm Local_spin_resolved_densities0.0001000000.txt
rm a.out
date
