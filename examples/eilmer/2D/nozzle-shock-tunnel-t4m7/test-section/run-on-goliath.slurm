#!/bin/bash
#SBATCH --job-name=t4m7-test-section
#SBATCH --nodes=1
#SBATCH --ntasks=22

module load mpi/openmpi-x86_64
mpirun e4mpi --job=t4m7-test-section --run  > LOGFILE_RUN

