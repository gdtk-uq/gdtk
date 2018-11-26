#!/bin/bash
prep-gas nitrogen-2sp.inp nitrogen-2sp.lua
prep-chem --compact nitrogen-2sp.lua nitrogen-2sp-2r.lua e4-chem.lua
prep-cuda-gpu-chem-kernel.py
e4shared --prep --job=n90
