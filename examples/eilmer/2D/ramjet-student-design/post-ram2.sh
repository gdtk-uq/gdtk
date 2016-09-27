#!/bin/bash
# post-ram2.sh
e4shared --post --job=ram2 --tindx-plot=all --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
