#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=exptube
mpirun -np 3 e4mpi --run --job=exptube --verbosity=1
e4shared --post --job=exptube --vtk-xml --tindx-plot=all \
	 --add-vars="mach,pitot,total-p,total-h"
