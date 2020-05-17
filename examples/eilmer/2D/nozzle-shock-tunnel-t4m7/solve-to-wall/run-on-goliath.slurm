#!/bin/bash
#SBATCH --job-name=t4m7
#SBATCH --nodes=1
#SBATCH --ntasks=4

module load mpi/openmpi-x86_64
mpirun e4mpi --job=t4m7 --run  > LOGFILE_RUN

