#!/bin/bash
e4shared --post --job=vtx --tindx-plot=last \
    --ref-soln=udf-vortex-flow.lua \
    --norms="p"

e4shared --post --job=vtx --tindx-plot=last \
    --ref-soln=udf-vortex-flow.lua \
    --norms="T"

