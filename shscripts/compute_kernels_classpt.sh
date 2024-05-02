#!/bin/bash
#SBATCH --job-name=compute_kernels
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=1G
#SBATCH --output=out/%x.%j.out
#SBATCH --error=err/%x.%j.err

cd /home/rneveux/bispectrum/theory/

module load anaconda

f=$(($2+124))
for i in `seq $2 $f`
do
python compute_kernels_classpt.py -redshift $1 -cosmo $i -directory $3 -Omfid $4
done