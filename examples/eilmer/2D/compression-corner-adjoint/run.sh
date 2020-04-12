#!/bin/bash

# prep
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=comp-corner

# run flow solver
e4-nk-shared --job=comp-corner

# post
e4shared --job=comp-corner --post --tindx-plot=last --vtk-xml

# parameterise ramp surface
#e4ssc --job=comp-corner --parameterise-surfaces --adjoint-verification=true

# run adjoint solver
e4ssc --job=comp-corner --adjoint-method --adjoint-verification=true
