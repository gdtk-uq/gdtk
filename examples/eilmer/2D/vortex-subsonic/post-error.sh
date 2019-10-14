#!/bin/bash
# post-error.sh
e4shared --post --job=vortex --tindx-plot=last \
         --ref-soln=vortex-flow-spec.lua \
         --norms="rho,p,T,vel.x,vel.y"

