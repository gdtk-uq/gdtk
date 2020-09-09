#!/bin/bash
# post.sh
e4shared --post --job=t4m7 --tindx-plot=last --vtk-xml \
	        --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=t4m7 --tindx-plot=last \
   --add-vars="mach,pitot,total-p,total-h" \
   --slice-list="112,18:19,:,0;113,18:19,:,0;114,18:19,:,0;115,18:19,:,0" \
   --output-file="extract-slice.data"
