#!/bin/bash
prep-gas gas.inp air-5sp-2T.lua
prep-chem air-5sp-2T.lua GuptaEtAl-air-2T.lua air-5sp-6r-2T.lua
e4shared --job=nonaka --prep
