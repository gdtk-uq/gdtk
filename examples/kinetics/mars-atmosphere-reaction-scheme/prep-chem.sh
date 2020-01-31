#!/bin/bash
# prep-chem.sh
prep-gas mars-atm-with-ions.inp mars-atm-with-ions.lua
prep-chem mars-atm-with-ions.lua Park_et_al_mars_atm_reactions.lua mars-atm-with-ions-chemistry.lua
