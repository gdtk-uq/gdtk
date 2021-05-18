#!/bin/bash
prep-gas air-5sp-2T.inp air-5sp-2T.gas
prep-chem air-5sp-2T.gas GuptaEtAl-air-2T.lua air-5sp-6r-2T.chem
prep-kinetics air-5sp-2T.gas air-5sp-6r-2T.chem air-VT-energy-exchange.lua air-VT.exch
e4shared --job=nonaka --prep
