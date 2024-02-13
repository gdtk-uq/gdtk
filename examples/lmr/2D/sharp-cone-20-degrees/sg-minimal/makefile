LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep: gas grid init-steady

gas: ideal-air.lua
	prep-gas ideal-air.lua ideal-air.gas

grid: grid.lua
	lmr prep-grid --job=grid.lua

init-steady: steady.lua
	lmr prep-sim --job=steady.lua

init-transient: transient.lua
	lmr prep-sim --job=transient.lua

run:
	lmr run

post: vtk

vtk:
	lmr snapshot2vtk --all

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
