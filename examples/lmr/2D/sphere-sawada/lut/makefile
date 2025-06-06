LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep:
	cp $(DGD_REPO)/src/gas/sample-data/air-5sp-1T-input.lua ./
	cp $(DGD_REPO)/src/gas/sample-data/air-5sp-eq.lua ./
	lmr prep-gas -i air-5sp-1T-input.lua -o air-5sp-1T.lua
	build-uniform-lut --gas-model=air-5sp-eq.lua --table-name=air5
	lmr prep-grid
	lmr prep-sim

run:
	mpirun -np 4 lmr-mpi-run

post:
	lmr snapshot2vtk --all

standoff:
	lmr custom-script --job=shock-standoff.lua

clean:
	- rm -rf $(LMR_DIR) lut-air5.lua air-5sp-1T-input.lua air-5sp-1T.lua air-5sp-eq.lua
