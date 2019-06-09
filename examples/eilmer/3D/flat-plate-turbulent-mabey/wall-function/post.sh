#!/bin/bash
# run.sh

# Extracts Paraview files for flowfield visualisation
#e4shared --post --job=mabey_3D --vtk-xml --tindx-plot=all \
# --add-vars="mach,pitot,total-p,total-h"

# Extracts slices at x=0.368m. This will allow us to
# examine the boundary layer profiles at these locations.
# Have looked in viscous_data.dat file to determine which 
# i-index to select.
e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_1.dat" \
  --slice-list="6,28,:,0;7,28,:,0" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_2.dat" \
  --slice-list="6,28,:,1;7,28,:,1" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_3.dat" \
  --slice-list="6,28,:,2;7,28,:,2" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_4.dat" \
  --slice-list="6,28,:,3;7,28,:,3" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_5.dat" \
  --slice-list="6,28,:,4;7,28,:,4" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_6.dat" \
  --slice-list="6,28,:,5;7,28,:,5" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_7.dat" \
  --slice-list="6,28,:,6;7,28,:,6" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_8.dat" \
  --slice-list="6,28,:,7;7,28,:,7" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_9.dat" \
  --slice-list="6,28,:,8;7,28,:,8" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_10.dat" \
  --slice-list="6,28,:,9;7,28,:,9" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-x-368mm_11.dat" \
  --slice-list="6,28,:,10;7,28,:,10" \
  --add-vars="mach,pitot,total-p,total-h"

# Extracts a slice of the nearest cells (to the wall) along 
# the plate. This will allow us to examine viscous properties,
# like skin friction and y+ values.
e4shared --post --job=mabey_3D --output-file="mabey-y-wall_1.dat" \
  --slice-list="1,:,42,0;3,:,42,0;5,:,42,0;7,:,42,0" \
  --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=mabey_3D --output-file="mabey-y-wall_2.dat" \
  --slice-list="1,:,42,1;3,:,42,1;5,:,42,1;7,:,42,1" \
  --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=mabey_3D --output-file="mabey-y-wall_3.dat" \
  --slice-list="1,:,42,2;3,:,42,2;5,:,42,2;7,:,42,2" \
  --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=mabey_3D --output-file="mabey-y-wall_4.dat" \
  --slice-list="1,:,42,3;3,:,42,3;5,:,42,3;7,:,42,3" \
  --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=mabey_3D --output-file="mabey-y-wall_5.dat" \
  --slice-list="1,:,42,4;3,:,42,4;5,:,42,4;7,:,42,4" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey_3D --output-file="mabey-y-wall_6.dat" \
  --slice-list="1,:,42,5;3,:,42,5;5,:,42,5;7,:,42,5" \
  --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=mabey_3D --output-file="mabey-y-wall_7.dat" \
  --slice-list="1,:,42,6;3,:,42,6;5,:,42,6;7,:,42,6" \
  --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=mabey_3D --output-file="mabey-y-wall_8.dat" \
  --slice-list="1,:,42,7;3,:,42,7;5,:,42,7;7,:,42,7" \
  --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=mabey_3D --output-file="mabey-y-wall_9.dat" \
  --slice-list="1,:,42,8;3,:,42,8;5,:,42,8;7,:,42,8" \
  --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=mabey_3D --output-file="mabey-y-wall_10.dat" \
  --slice-list="1,:,42,9;3,:,42,9;5,:,42,9;7,:,42,9" \
  --add-vars="mach,pitot,total-p,total-h"
e4shared --post --job=mabey_3D --output-file="mabey-y-wall_11.dat" \
  --slice-list="1,:,42,10;3,:,42,10;5,:,42,10;7,:,42,10" \
  --add-vars="mach,pitot,total-p,total-h"


# Runs python script to compute viscous parameters.
#echo "Computing viscous parameters .."
#python compute_viscous_data_along_wall_wallFunction.py

#Perform van Driest transformation to obtain dimensionless parameters.
#echo "Perform van Driest transformation .."
#python compute_viscous_data_along_wall_wallFunction.py

######################## Transfer plotting files to combined folder.
#echo "Transfering .dat files.."
#./transfer.sh

########################  Plots numerical results against experimental data.
#echo "Plotting results .."
#./plot.sh

echo "Postprocessing completed."
