LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run-steady"
	@echo -e "   post"

prep: gas grid init

gas:
	echo "model = 'ideal gas'; species = {'air'}" > ideal-air.lua 
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas

grid:
	cd ../../larc-grids; lmr prep-grid --job=lmr-coarsest-grid.lua
	- mkdir -p lmrsim/grid
	rsync -a --delete ../../larc-grids/lmrsim/grid/ lmrsim/grid/
	- rm -rf ../../larc-grids/lmrsim
	lmr structured2unstructured

init: job.lua
	lmr prep-sim

run:
	mpirun -n 10 lmrZ-mpi-run

post: vtk

vtk:
	lmr snapshot2vtk --all

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.lua ideal-air.gas
