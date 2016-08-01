#!/bin/bash

prep-gas Stanford-2011-species.inp Stanford-2011-gas-model.lua
prep-chem Stanford-2011-gas-model.lua Stanford-2011.lua Stanford-2011-reac-file.lua

sed -i '1 s/.*/spFile = "Stanford-2011-gas-model.lua"/' ignition-delay.lua
sed -i '2 s/.*/reacFile = "Stanford-2011-reac-file.lua"/' ignition-delay.lua
sed -i '3 s/.*/outFile = "Stanford-ignition-delay.dat"/' ignition-delay.lua

gas-calc ignition-delay.lua

