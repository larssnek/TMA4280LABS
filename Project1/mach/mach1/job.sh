#!/bin/bash

#PBS -N mach1
#PBS -q training
#PBS -A imf_lille-tma4280
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:15:00
#PBS -l select=2:ncpus=20:mpiprocs=16

cd $PBS_O_WORKDIR
module load gcc
module load openmpi
mpiexec ./mach1
