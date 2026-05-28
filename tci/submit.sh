#!/bin/sh

#SBATCH -t 48:00:00
#SBATCH -o sbatch.out
#SBATCH -e sbatch.err
#SBATCH -c 10
#SBATCH -N 1
#SBATCH --mem=60000
#SBATCH --account=berkelbach

# export lib path for user specific lib

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

python tci_mol.py > stdout
