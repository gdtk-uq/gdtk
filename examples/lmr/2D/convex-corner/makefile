LMR_DIRS := config grid snapshots vtk diagnostics

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

links:
	- ln -s lmr-job.lua lmr-grid.lua
	- ln -s lmr-job.lua lmr-flow.lua

prep: links prep-gas grid init

prep-gas: gm-ideal-air.inp
	prep-gas gm-ideal-air.inp ideal-air.gas

.PHONY: grid
grid:
	lmr prep-grid

init:
	lmr prep-flow --mode=steady

run:
	lmr run-steady

post: vtk

.PHONY: vtk
vtk:
	lmr snapshot2vtk --all

clean:
	rm -rf $(LMR_DIRS)

