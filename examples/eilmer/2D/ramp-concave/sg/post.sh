#!/bin/bash
# post.sh
e4shared --post --job=cubic-ramp --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
