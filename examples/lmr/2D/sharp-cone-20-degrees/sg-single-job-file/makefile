LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep: gas grid init

gas: ideal-air.lua
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas

grid: job.lua
	lmr prep-grid

init: job.lua
	lmr prep-sim

run:
	lmr run

post: vtk

vtk:
	lmr snapshot2vtk --all

test:
	PYTHONDONTWRITEBYTEFILE=1 pytest -v

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
