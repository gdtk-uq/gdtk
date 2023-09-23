#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=ns-mms-2D

source $HOME/cupy/bin/activate
make verify

