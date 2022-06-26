#!/bin/bash

# Copy Gupta et al across, rather than store yet another copy.
cp ../air-chemistry-1T/GuptaEtAl-air-reactions.lua .

# Ensure that 5-species for air is selected
sed -i "s/^SUBSET_SELECTION.*/SUBSET_SELECTION = '5-species-air'/" GuptaEtAl-air-reactions.lua

# Then prep gas and chemistry in the usual manner
prep-gas air-5sp-He-1T.lua air-5sp-He-1T.gas
prep-chem air-5sp-He-1T.gas GuptaEtAl-air-reactions.lua GuptaEtAl-air-reactions.chem


