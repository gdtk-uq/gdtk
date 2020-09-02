#!/bin/bash
# prep.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=annulus
e4shared --post --job=annulus --vtk-xml --tindx-plot=0 \
	 --add-vars="mach,pitot,total-p,total-h,nrf"
