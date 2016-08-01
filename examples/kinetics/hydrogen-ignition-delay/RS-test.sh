#!/bin/bash

prep-gas Rogers-Schexnayder-species.inp Rogers-Schexnayder-gas-model.lua
prep-chem Rogers-Schexnayder-gas-model.lua Rogers-Schexnayder.lua Rogers-Schexnayder-reac-file.lua

sed -i '1 s/.*/spFile = "Rogers-Schexnayder-gas-model.lua"/' ignition-delay.lua
sed -i '2 s/.*/reacFile = "Rogers-Schexnayder-reac-file.lua"/' ignition-delay.lua
sed -i '3 s/.*/outFile = "RS-ignition-delay.dat"/' ignition-delay.lua

gas-calc ignition-delay.lua

