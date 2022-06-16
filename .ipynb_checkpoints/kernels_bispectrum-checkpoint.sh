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
python kernels_bispectrum.py -cosmo $i  -redshift .8 -ells $l -directory /home/rneveux/kernels_EFT/bk/z0.8/Omfid.31377/kernels/ -AP True -k_list /home/rneveux/kernels_EFT/bk/z0.8/Omfid.31377/k_data.txt
done
done