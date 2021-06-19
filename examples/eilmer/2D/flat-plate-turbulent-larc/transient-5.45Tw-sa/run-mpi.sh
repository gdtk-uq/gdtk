#!/bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=larc
mpirun -np 6 --oversubscribe e4mpi --job=larc --run --verbosity=1
e4shared --post --job=larc --vtk-xml --add-vars="mach,pitot,total-p,total-h" --tindx-plot=last
python3 compute_drag.py larc
