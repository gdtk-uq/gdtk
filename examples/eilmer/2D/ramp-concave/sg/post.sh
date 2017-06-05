#!/bin/bash
# post.sh
e4shared --post --job=cubic-ramp --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"

gnuplot surface-pressure.gnuplot

awk -f scale-heat-flux.awk ./loads/t10-loads.dat > stanton.data
gnuplot surface-heat-transfer.gnuplot
