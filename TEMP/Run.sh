#!/bin/bash
#SBATCH -p all-nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60000mb
#SBATCH --job-name=HF_LL
module load gcc/9.2
module load mkl
module load ld-path
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

hostname
date
time ./k_space_SelfConsistency Generic2dLattice input_run.inp > out_run.txt
date
