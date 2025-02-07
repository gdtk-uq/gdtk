# prep-gas-models.sh
# PJ 2025-02-06

cp ${DGD_REPO}/src/gas/sample-data/air-5sp-1T-input.lua .
cp ${DGD_REPO}/src/gas/sample-data/air-5sp-eq.lua .
prep-gas air-5sp-1T-input.lua air-5sp-1T.lua

cp ${DGD_REPO}/src/gas/sample-data/co2-5sp-1T-input.lua .
cp ${DGD_REPO}/src/gas/sample-data/co2-5sp-eq.lua .
prep-gas co2-5sp-1T-input.lua co2-5sp-1T.lua

cp ${DGD_REPO}/src/gas/sample-data/cea-co2-gas-model.lua .
