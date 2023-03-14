#!/bin/bash
# run_partiioner.sh

rm block*
e4shared --custom-script --script-file=su2-grid-gen.lua
ugrid_partition ramp15.su2 mapped_cells 8 2
mv mapped_cells ../
