LMR_DIR := lmrsim

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep:
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-5sp-gas-model.lua .
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/GuptaEtAl-air-reactions-2T.lua .
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-energy-exchange.lua .
	prep-gas air-5sp-gas-model.lua air-5sp-2T.gas
ifeq (${PLATFORM}, macosx)
	# instruct mac to make no backup file, that's the empty string up front
	sed -i '' '/^SUBSET_SELECTION =/s/11/5/' GuptaEtAl-air-reactions-2T.lua
else
	sed -i '/^SUBSET_SELECTION =/s/11/5/' GuptaEtAl-air-reactions-2T.lua
endif
	prep-chem air-5sp-2T.gas GuptaEtAl-air-reactions-2T.lua air-5sp-6r-2T.chem
	prep-kinetics air-5sp-2T.gas air-5sp-6r-2T.chem air-energy-exchange.lua air-VT.exch
	lmr prep-grid --job=grid.lua

prep-transient: transient.lua
	lmr prep-sim --job=transient.lua

prep-steady: steady.lua
	lmr prep-sim --job=steady.lua

run:
	lmr run

post:
	lmr snapshot2vtk --all

clean:
	- rm -rf $(LMR_DIR)
	- rm air-5sp-gas-model.lua air-5sp-2T.gas GuptaEtAl-air-reactions-2T.lua
	- rm air-energy-exchange.lua air-5sp-6r-2T.chem air-VT.exch
