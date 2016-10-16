#!/bin/bash
# prep.sh
prep-gas nitrogen-2sp.inp nitrogen-2sp.lua
prep-chem nitrogen-2sp.lua nitrogen-2sp-2r.lua e4-chem.lua
e4shared --prep --job=cyl
