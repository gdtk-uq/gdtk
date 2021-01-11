#!/bin/bash
prep-gas N2-2T-diss.inp nitrogen-gas-model.lua
prep-chem nitrogen-gas-model.lua nitrogen-2sp-2T-2r.lua 2T-dissociating-N2.lua
e4shared --job=n90 --prep
