#!/bin/bash

# Copy Gupta et al across, rather than store yet another copy.
cp ../air-chemistry-2T/GuptaEtAl-air-reactions-2T.lua .

# Ensure that 5-species for air is selected
sed -i "s/^SUBSET_SELECTION.*/SUBSET_SELECTION = '5-species-air'/" GuptaEtAl-air-reactions-2T.lua

# Copy over energy exchange file
cp ../air-chemistry-2T/air-energy-exchange.lua .

# Ensure that electrons are not included 
sed -i "s/^WITH_ELECTRONS.*/WITH_ELECTRONS = false/" air-energy-exchange.lua

# Then prep gas, chemistry and energy exchange in the usual manner
prep-gas air-5sp-He-2T.lua air-5sp-He-2T.gas
prep-chem air-5sp-He-2T.gas GuptaEtAl-air-reactions-2T.lua GuptaEtAl-air-reactions-2T.chem
prep-kinetics air-5sp-He-2T.gas air-energy-exchange.lua air-energy-exchange.kin

