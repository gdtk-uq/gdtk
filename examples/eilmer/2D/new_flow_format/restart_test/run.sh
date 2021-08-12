#!/bin/bash
# run.sh

prep-gas ideal-air.inp ideal-air-gas-model.lua

e4shared --prep --job=ffs
e4shared --run --job=ffs --verbosity=1 --max-cpus=3
e4shared --post --job=ffs --tindx-plot=all --vtk-xml

e4shared --prep --job=ffs_restart_1
e4shared --run --job=ffs_restart_1 --verbosity=1 --max-cpus=3
e4shared --post --job=ffs_restart_1 --tindx-plot=all --vtk-xml

e4shared --prep --job=ffs_restart_2
e4shared --run --job=ffs_restart_2 --verbosity=1 --max-cpus=3
e4shared --post --job=ffs_restart_2 --tindx-plot=all --vtk-xml

