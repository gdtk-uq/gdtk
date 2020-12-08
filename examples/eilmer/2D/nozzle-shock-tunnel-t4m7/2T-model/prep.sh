#!/bin/bash
# prep.sh
# Do the preparation of grid and initial flow state, only.
# Note that we use the debug flavour of the code for the CEA gas model
# that gets used in the calculation of the throat conditions.
e4shared --prep --job=t4m7_noneq

# Postprocessing uses only the LUT gas model.
e4shared --post --job=t4m7_noneq --tindx-plot=0 --vtk-xml

