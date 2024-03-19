LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep: ideal-air.lua job.lua
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas
	lmr prep-grid
	lmr prep-sim

run:
	lmr run

post:
	lmr snapshot2vtk --all
	lmr slice-flow --slice-list="0,:,0,0;2,:,0,0;4,:,0,0" --output=stag-line.data

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
	- rm stag-line.data
