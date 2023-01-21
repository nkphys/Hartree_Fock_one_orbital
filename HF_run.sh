#!/bin/bash
#SBATCH --job-name=HFJob
#SBATCH --account ACF-UTK0038
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=condo-dagotto
#SBATCH --qos=condo
#SBATCH --time=100:00:00
#SBATCH -o ParamsJob.%j


hostname
module swap PE-intel PE-gnu
module load lapack
module load hdf5
module load openBLAS
module load mkl
date
time  ./k_space_SelfConsistency HoneycombLattice_MO input.inp > out_run.txt
date
