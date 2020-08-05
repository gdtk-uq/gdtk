#!/bin/bash
# run.sh
e4shared --post --job=t4m4b --tindx-plot=last --vtk-xml \
   --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=t4m4b --tindx-plot=last \
   --add-vars="mach,pitot,total-p,total-h" \
   --slice-list="60,$,:,0;61,$,:,0" --output-file="nozzle-exit-eilmer4.data"
