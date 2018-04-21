#!/bin/bash
# prep-chem.sh
prep-gas nitrogen-2sp.inp nitrogen-2sp.lua
prep-chem nitrogen-2sp.lua nitrogen-2sp-2r.lua nitrogen-chemistry.lua

prep-gas air-5sp.inp air-5sp.lua
prep-chem air-5sp.lua GuptaEtAl-air-reactions.lua air-chemistry.lua
