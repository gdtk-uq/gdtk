#!/bin/bash
# run.sh
# time e4shared --run --job=swbli --verbosity=1 --max-cpus=4 --report-residuals
time mpirun -np 4 e4mpi --run --job=swbli --verbosity=1 --report-residuals
