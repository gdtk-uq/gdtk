#!/bin/bash
# run-fine-restart.sh
e4shared --prep --job=cone20-fine-restart
e4shared --run --job=cone20-fine-restart --verbosity=1 --max-cpus=2
e4shared --post --job=cone20-fine-restart --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
