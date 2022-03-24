#! /bin/bash
# run.sh
e4shared --job=sphere-cone --run
e4shared --job=sphere-cone --post --vtk-xml --tindx-plot=all
