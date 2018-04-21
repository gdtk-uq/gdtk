#!/bin/bash

prep-gas drm19-prep-gas.inp drm19-gas-model.lua
prep-chem drm19-gas-model.lua drm19-prep-chem.inp drm19-reaction-scheme.lua
gas-calc fixed-volume-reactor.lua > output.data
gnuplot plot-comparison.gplot
ps2pdf -dEPSCrop drm19-fixed-volume-reactor.eps

