#!/bin/bash
prep-gas O2.inp O2-gas-model.lua
gas-calc thermo-curves-for-O2.lua
gnuplot plot-data.gplot
ps2pdf -dEPSCrop O2-thermo.eps

