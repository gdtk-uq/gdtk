#!/bin/bash
# run.sh

prep-gas ideal-N2.inp ideal-N2-gas-model.lua
e4shared --prep --job=mallinson
e4shared --run --job=mallinson --verbosity=1 --max-cpus=4
