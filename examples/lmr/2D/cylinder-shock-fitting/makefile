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

post: vtk slice

vtk:
	lmr snapshot2vtk --all

slice:
	lmr slice-flow --slice-list="0,:,0,0;2,:,0,0;4,:,0,0"

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
