#!/bin/bash
#SBATCH --job-name=kernel_computation_bk
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --output=out/%x.%j.out
#SBATCH --error=err/%x.%j.err

cd /home/rneveux/bispectrum/theory/

module load anaconda

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/rneveux/cuba

for j in `seq 0 51`
do
python compute_kernels.py -spectrum_part $2 -kernel $j -ells $1 -estimator bk -redshift .8 -directory /home/rneveux/kernels/bk/z.8/k_1_20/no_f/ -f_fit oui
done