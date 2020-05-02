#!/bin/bash
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/air-5sp-1T.inp ./air-5sp.inp
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua .
prep-gas air-5sp.inp air-5sp.lua
prep-chem air-5sp.lua GuptaEtAl-air-reactions.lua air-5sp-6r.lua
