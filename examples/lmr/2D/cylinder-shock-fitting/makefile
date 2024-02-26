LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep: gas grid init-transient

gas: ideal-air.lua
	prep-gas ideal-air.lua ideal-air.gas

grid: job.lua
	lmr prep-grid

init-transient: job.lua
	lmr prep-sim

run:
	lmr run

post: vtk

vtk:
	lmr snapshot2vtk --all

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
