#!/bin/bash
#PBS -l walltime=00:15:00
#PBS -l nodes=1:ppn=20
#PBS -W group_list=newriver
#PBS -q normal_q
#PBS -A CMDA3634
#Change to the directory from which the job was submitted
cd $PBS_O_WORKDIR
#Load modules
module purge; module load gcc openmpi
#May not be necessary if the program is already built
gcc -O3 -o mandelbrot mandelbrot.c -lm -fopenmp 
#Run
#echo "Run mpiBTNPOT with $PBS_NP cores"
for i in `seq 1 20`
do
./mandelbrot 4096 4096 $i
done
exit;
