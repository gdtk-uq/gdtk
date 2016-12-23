#!/bin/bash
prep-gas thermally-perfect-N2-O2.inp thermally-perfect-N2-O2.lua
gas-calc transport-properties-for-air.lua
gnuplot plot-data.gplot
ps2pdf -dEPSCrop air-trans.eps

