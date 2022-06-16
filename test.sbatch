#!/bin/bash
#SBATCH --ntasks=10
#SBATCH --time=10:00:00
#SBATCH --mem=1G

for i in `seq $1 10 $(($1 + 90))`
do
echo 'ok'
echo $i
#for j in `seq 0 51`
#do
mpirun python test.py $i $j
#done
done
