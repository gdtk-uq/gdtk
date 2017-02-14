#!/bin/bash
# runB.sh
e4shared --post --job=cone20 --extract-line="0,0,0,0,1,0,100" --output-file="profile.data"
e4shared --prep --job=cone20B
e4shared --run --job=cone20B --verbosity=1 --max-cpus=2
e4shared --post --job=cone20B --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
