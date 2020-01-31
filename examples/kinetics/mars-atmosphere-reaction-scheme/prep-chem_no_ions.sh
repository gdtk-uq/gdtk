#!/bin/bash
# prep-chem.sh
prep-gas mars-atm.inp mars-atm.lua
prep-chem mars-atm.lua Park_et_al_mars_atm_reactions.lua mars-atm-chemistry.lua
