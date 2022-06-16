#!/bin/bash
#SBATCH --job-name=compute_bk
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --mem=1G
#SBATCH --output=out/%x.%j.out
#SBATCH --error=err/%x.%j.err

cd /home/rneveux/bispectrum/theory/

module load anaconda

python computation_bk.py -ells $1