LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run-steady"
	@echo -e "   post"

prep: gas grid init

gas: ideal-air.lua
	prep-gas ideal-air.lua ideal-air.gas

grid: grid.lua
	lmr prep-grid --job=grid.lua

init: steady.lua
	lmr prep-sim --job=steady.lua

run-steady:
	lmr run

post: vtk

vtk:
	lmr snapshot2vtk --all

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
