# run.sh
# Expansion tube demo.
#
prep-gas ideal-air.inp ideal-air-gas-model.lua
l1d4-prep --job=exptube
l1d4 --run-simulation --job=exptube
