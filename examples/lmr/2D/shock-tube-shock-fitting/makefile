LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep: ideal-air.lua job.lua
	prep-gas ideal-air.lua ideal-air.gas
	lmr prep-grid
	lmr prep-sim

run:
	lmr run

post:
	lmr snapshot2vtk --all

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
