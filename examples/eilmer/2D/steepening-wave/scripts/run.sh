#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=swp
e4shared --run --job=swp --verbosity=1
