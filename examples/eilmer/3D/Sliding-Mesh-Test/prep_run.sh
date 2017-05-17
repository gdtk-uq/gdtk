#! /bin/sh
# run_Sliding_Mesh_test.sh
# Test Case for Developing a sliding mesh interface

# Create gas model
prep-gas ideal-air.inp ideal-air-gas-model.lua

# Create the grid.
e4shared --prep --job=Sliding_Mesh_Test

# run the code
#sleep 10 && gnuplot live-residuals.gplot &
e4shared --run --job=Sliding_Mesh_Test --verbosity=1 --max-cpus=1
