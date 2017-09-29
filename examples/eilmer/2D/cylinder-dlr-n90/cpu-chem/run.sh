#!/bin/bash
# run.sh
e4shared --prep --job=n90
e4shared --run --job=n90 --verbosity=1 --max-cpus=4
e4shared --post --job=n90 --tindx-plot=all --vtk-xml
