#!/bin/bash
#SBATCH --job-name=n90
#SBATCH --gres=gpu:1
#SBATCH --mincpus=1
#SBATCH --partition=gpu
#SBATCH --exclusive
#SBATCH --nodelist="c4130-2"

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. eilmer-cuda-gpu-chem --run --job=n90 --max-cpus=1 > LOGFILE 2> ERRFILE


