#!/bin/bash
# run.sh
e4shared --run --job=annulus --verbosity=1
e4shared --post --job=annulus --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h,nrf"
