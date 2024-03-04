LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep:
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-5sp-gas-model.lua .
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/GuptaEtAl-air-reactions-2T.lua .
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-energy-exchange.lua .
	sed -i '/^SUBSET_SELECTION =/s/11/5/' GuptaEtAl-air-reactions-2T.lua
	prep-gas air-5sp-gas-model.lua air-5sp-2T.gas
	prep-chem air-5sp-2T.gas GuptaEtAl-air-reactions-2T.lua air-5sp-6r-2T.chem
	prep-kinetics air-5sp-2T.gas air-5sp-6r-2T.chem air-energy-exchange.lua air-VT.exch
	lmr prep-grid
	lmr prep-sim

run:
	lmr run

post:
	lmr snapshot2vtk --all

clean:
	- rm -rf $(LMR_DIR)
	- rm air-5sp-gas-model.lua air-5sp-2T.gas GuptaEtAl-air-reactions-2T.lua
	- rm air-energy-exchange.lua air-5sp-6r-2T.chem air-VT.exch
