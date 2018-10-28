#!/bin/bash
# post.sh

job=eklund

# xh_-10 = -3.175e-3 m
# xh_-01 = -0.3175e-3 m
# xh_17 = 5.3975e-3 m
# xh_30 = 9.525e-3 m
# xh_39 = 12.3825e-3 m
# xh_67 = 21.2725e-3 m
# xh_108 = 34.29e-3 m

# Extracts Paraview files for flowfield visualisation
e4shared --post --job=$job --vtk-xml --add-vars="mach,pitot,total-p,total-h"

# Boundary layer slices to compare with experimental data
e4shared --post --job=$job --output-file="xh_-10.dat" \
  --slice-list="15,18,:,:;16,18,:,:" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=$job --output-file="xh_-01.dat" \
  --slice-list="15,18,:,:;16,18,:,:" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=$job --output-file="xh_17.dat" \
  --slice-list="1,21,:,:;7,21,:,:;8,21,:,:" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=$job --output-file="xh_30.dat" \
  --slice-list="2,14,:,:;9,14,:,:;10,14,:,:" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=$job --output-file="xh_39.dat" \
  --slice-list="2,25,:,:;9,25,:,:;10,25,:,:" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=$job --output-file="xh_67.dat" \
  --slice-list="3,24,:,:;11,24,:,:;12,24,:,:" \
  --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=$job --output-file="xh_108.dat" \
  --slice-list="4,$,:,:;13,$,:,:;14,$,:,:" \
  --add-vars="mach,pitot,total-p,total-h"

echo "Postprocessing completed."
