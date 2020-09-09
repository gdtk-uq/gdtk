#!/bin/bash
# post.sh
e4shared --post --job=t4m7-test-section --tindx-plot=all --vtk-xml \
	       --add-vars="mach,pitot,total-p,total-h"

