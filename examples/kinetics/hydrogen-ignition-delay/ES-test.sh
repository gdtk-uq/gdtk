#!/bin/bash

prep-gas Evans-Schexnayder-species.inp Evans-Schexnayder-gas-model.lua
prep-chem Evans-Schexnayder-gas-model.lua Evans-Schexnayder.lua Evans-Schexnayder-reac-file.lua

sed -i '1 s/.*/spFile = "Evans-Schexnayder-gas-model.lua"/' ignition-delay.lua
sed -i '2 s/.*/reacFile = "Evans-Schexnayder-reac-file.lua"/' ignition-delay.lua
sed -i '3 s/.*/outFile = "ES-ignition-delay.dat"/' ignition-delay.lua

gas-calc ignition-delay.lua

