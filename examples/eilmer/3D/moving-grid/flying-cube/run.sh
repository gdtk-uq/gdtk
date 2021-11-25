#! /bin/bash
# run.sh
e4shared --job=cube --run
e4shared --job=cube --post --vtk-xml --tindx-plot=all
