#!/bin/bash
# run.sh
e4shared --prep --job=cyl-sf
e4shared --run --job=cyl-sf --verbosity=1 --max-cpus=4
e4shared --post --job=cyl-sf --tindx-plot=all --vtk-xml
