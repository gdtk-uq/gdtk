#!/bin/bash
# Extract the stagnation line data from the steady flow field.
e4shared --post --job=n90 --output-file=n90-stag-prof-KH.data --tindx-plot=5 \
    --slice-list="0,:,1,0"


