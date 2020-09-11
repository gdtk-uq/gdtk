#!/bin/bash

# I've manually modified thermally-perfect-N2-O2-mix to have the
# same lewis numbers for N2 and O2, so we can compare to the 
# analytical solution. So don't regenerate it if that's 
# what you're trying to do.
#prep-gas gas-model.inp thermally-perfect-N2-O2-mix.lua
e4shared --job=diffusion --prep

