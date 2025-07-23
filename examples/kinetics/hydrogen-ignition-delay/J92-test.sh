#!/bin/bash

prep-gas Jachimowski-1992-species.inp Jachimowski-1992-gas-model.lua
prep-chem Jachimowski-1992-gas-model.lua Jachimowski-1992.lua Jachimowski-1992-reac-file.lua

sed -i '1 s/.*/spFile = "Jachimowski-1992-gas-model.lua"/' ignition-delay.lua
sed -i '2 s/.*/reacFile = "Jachimowski-1992-reac-file.lua"/' ignition-delay.lua
sed -i '3 s/.*/outFile = "J92-ignition-delay.dat"/' ignition-delay.lua

gas-calc ignition-delay.lua

