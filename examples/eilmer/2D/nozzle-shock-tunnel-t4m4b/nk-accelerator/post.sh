#!/bin/bash
# run.sh
e4shared --post --job=t4m4b --tindx-plot=last --vtk-xml \
   --add-vars="mach,pitot,total-p,total-h"
