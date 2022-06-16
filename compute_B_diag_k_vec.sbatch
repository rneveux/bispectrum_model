#!/bin/bash
#SBATCH --job-name=compute_B_k_vec
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=1G
#SBATCH --output=out/%x.%j.out
#SBATCH --error=err/%x.%j.err

cd /home/rneveux/bispectrum/theory/

module load anaconda

python compute_B_diag_k_vec.py