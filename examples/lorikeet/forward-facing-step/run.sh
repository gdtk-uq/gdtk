#! /usr/bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
lrkt-prep --job=ffs
lrkt-run --job=ffs --max-cpus=7
lrkt-post --job=ffs --tindx=all --vtk-xml
