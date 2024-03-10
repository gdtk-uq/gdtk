LMR_DIR := lmrsim

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep:
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-5sp-gas-model.lua .
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/GuptaEtAl-air-reactions-2T.lua .
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-energy-exchange.lua .
	lmr prep-gas -i air-5sp-gas-model.lua -o air-5sp-2T.gas
ifeq (${PLATFORM}, macosx)
	# instruct mac to make no backup file, that's the empty string up front
	sed -i '' '/^SUBSET_SELECTION =/s/11/5/' GuptaEtAl-air-reactions-2T.lua
else
	sed -i '/^SUBSET_SELECTION =/s/11/5/' GuptaEtAl-air-reactions-2T.lua
endif
	lmr prep-reactions -g air-5sp-2T.gas -i GuptaEtAl-air-reactions-2T.lua -o air-5sp-6r-2T.chem
	lmr prep-energy-exchange -g air-5sp-2T.gas -r air-5sp-6r-2T.chem -i air-energy-exchange.lua -o air-VT.exch
	lmr prep-grid --job=grid.lua

prep-transient: transient.lua
	lmr prep-sim --job=transient.lua

prep-steady: steady.lua
	lmr prep-sim --job=steady.lua

run-transient:
	lmr run

run-steady:
	mpirun --use-hwthread-cpus -np 8 lmrZ-mpi-run

post:
	lmr snapshot2vtk --all

clean:
	- rm -rf $(LMR_DIR)
	- rm air-5sp-gas-model.lua air-5sp-2T.gas GuptaEtAl-air-reactions-2T.lua
	- rm air-energy-exchange.lua air-5sp-6r-2T.chem air-VT.exch
