#!/bin/bash
# run.sh
e4shared --run --job=t4m4 --verbosity=1 --max-cpus=4
e4shared --post --job=t4m4 --tindx-plot=last --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
