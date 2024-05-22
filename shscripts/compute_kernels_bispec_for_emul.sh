#!/bin/bash
#SBATCH --job-name=kernels_bispectrum
#SBATCH --ntasks=1
#SBATCH --time=50:00:00
#SBATCH --mem=1G
#SBATCH --output=out/%x.%j.out
#SBATCH --error=err/%x.%j.err

cd /home/rneveux/bispectrum/theory/

module load anaconda

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/rneveux/cuba

for l in 000 202
do
for i in `seq $1 $(($1 + 124))`
do
python pyscripts/compute_kernels_bispec_for_emul.py -cosmo $i  -redshift .61 -ells $l -directory /home/rneveux/kernels_EFT/bk/z0.61/Omfid.317/full_cosmo/kernels/ -Omfid .317 -k_list /home/rneveux/kernels_EFT/bk/z0.61/Omfid.317/k_data.txt
done
done