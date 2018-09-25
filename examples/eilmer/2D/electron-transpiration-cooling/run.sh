#!/bin/bash
# run.sh
# Script for alkandry2014 simulation. Preps and runs the
# simulation. Post-processing is performed by another script.
#
# Zachary J. Denman
# 2017-12-18
echo "========================================================================"
echo "================================ BEGIN ================================="
echo "========================================================================"

JOB_NAME="alkandry2014"

echo "Clearing old job files..."
make clean

# Generate gas model
echo "========================================================================"
echo "Generating gas model..."
echo "========================================================================"
prep-gas thermally-perfect-air-11sp.inp thermally-perfect-air-11sp.lua

# Run prep
echo "========================================================================"
echo "Preparing simulation..."
echo "========================================================================"
e4shared --prep --job=$JOB_NAME

# Run simulation
echo "========================================================================"
echo "Running simulation..."
echo "========================================================================"
e4shared --run --job=$JOB_NAME --verbosity=1 --max-cpus=8

# Post-processing
echo "========================================================================"
echo "Running post-processing..."
echo "========================================================================"
e4shared --post --job=$JOB_NAME --add-vars="mach" --tindx-plot=all --vtk-xml

echo "========================================================================"
echo "=============================== FINISHED ==============================="
echo "========================================================================"
