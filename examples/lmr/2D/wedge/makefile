LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep: gas grid init

gas: ideal-air.lua
	prep-gas ideal-air.lua ideal-air.gas

grid: job.lua
	lmr prep-grid

init:
	lmr prep-flow

run:
	lmr run

post: vtk

vtk:
	lmr snapshot2vtk --all

test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v
clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas

