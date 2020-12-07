#!/bin/bash
# post.sh
e4shared --post --job=t4m7_chem --tindx-plot=last --vtk-xml \
	       --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=t4m7_chem --tindx-plot=last \
   --add-vars="mach,pitot,total-p,total-h" \
   --slice-list="120,$,:,0;121,$,:,0;122,$,:,0;123,$,:,0" \
   --output-file="nozzle-exit-chem.data"


