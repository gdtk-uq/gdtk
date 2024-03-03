LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep: ideal-air.lua job.lua
	prep-gas ideal-air.lua ideal-air.gas
	lmr prep-grid
	lmr prep-sim

run:
	lmr run

post:
	lmr snapshot2vtk --all
	lmr slice-flow --output=sod_new.dat --slice-list="0,:,0,0"
	./plot.sh

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
