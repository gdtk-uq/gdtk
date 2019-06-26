#!/bin/bash
# run.sh
e4shared --run --job=icon --verbosity=1 --max-cpus=4
e4shared --post --job=icon --tindx-plot=all --vtk-xml \
         --add-vars="mach,pitot,total-p,total-h"

