#! /bin/bash
# run.sh
e4shared --job=sphere --run
e4shared --job=sphere --post --vtk-xml --tindx-plot=all
