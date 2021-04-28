#!/bin/bash
#SBATCH --job-name=cond1-grid2-inv
#SBATCH --nodes=1
#SBATCH --ntasks=1

l1d4 --run-simulation --job=x2_condition1 > LOGFILE
