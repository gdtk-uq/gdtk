#!/bin/bash
# run.sh
prep-gas nitrogen-2sp.inp nitrogen-2sp.lua
prep-chem --compact nitrogen-2sp.lua nitrogen-2sp-2r.lua e4-chem.lua
prep-opencl-gpu-chem-kernel.py
e4shared --prep --job=n90
eilmer-opencl-gpu-chem --run --job=n90 --verbosity=1 --max-cpus=1
e4shared --post --job=n90 --tindx-plot=all --vtk-xml
