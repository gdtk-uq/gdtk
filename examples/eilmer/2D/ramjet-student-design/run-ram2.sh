#!/bin/bash
# run-ram2.sh
e4shared --run --job=ram2 --verbosity=1 --max-cpus=4
e4shared --post --job=ram2 --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
