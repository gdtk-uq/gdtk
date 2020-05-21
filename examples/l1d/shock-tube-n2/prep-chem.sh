#!/bin/bash
# prep-chem.sh
prep-gas ideal-helium.inp ideal-helium-gas-model.lua
prep-gas nitrogen-2sp.inp nitrogen-2sp.lua
prep-chem nitrogen-2sp.lua nitrogen-2sp-2r-Keq.lua chem.lua
