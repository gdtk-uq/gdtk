LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep:
	cp ${DGD_REPO}/src/gas/sample-data/air-5sp-eq.lua ./
	cp ${DGD_REPO}/src/gas/sample-data/air-5sp-1T-input.lua ./
	lmr prep-gas -i air-5sp-1T-input.lua -o air-5sp-1T.lua
	lmr prep-grid
	lmr prep-sim

run:
	mpirun -np 4 lmr-mpi-run | tee log.txt

post:
	lmr snapshot2vtk --all

standoff:
	lmr custom-script --job=shock-standoff.lua

clean:
	- rm -rf $(LMR_DIR) air-5sp-1T-input.lua air-5sp-1T.lua air-5sp-eq.lua log.txt
