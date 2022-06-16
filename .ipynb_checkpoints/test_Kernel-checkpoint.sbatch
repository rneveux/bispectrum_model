#!/bin/bash
#SBATCH --job-name=test_Kernel
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
for j in `seq 0 51`
do
for k in tree SN PNG
do
python test_Kernel.py -cosmo $i -spectrum_part $k -kernel $j -ells $l -estimator bk -redshift .8 -directory /home/rneveux/kernels_EFT/bk/z0.8/lnAs/Omfid.31377/ -AP True
done
done
done
done
