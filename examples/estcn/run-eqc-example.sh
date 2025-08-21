#! /bin/bash
# Preapre and run an estcn example using the eqc equilibrium chemistry solver
cp ${DGD_REPO}/src/gas/sample-data/air-5sp-eq.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/air-5sp-1T.inp .
prep-gas air-5sp-1T.inp air-5sp-1T.lua
estcn --task=stn --gas=air-5sp-eq.lua --T1=300 --p1=270.0e3 --Vs=1313.03 --pe=18.858e6 --ar=157.8
