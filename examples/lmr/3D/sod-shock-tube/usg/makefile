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
	lmr extract-line --output=sod_new.dat --line-list="0.0,0.0,0.0,1.0,0.0,0.0,100"
	./plot.sh

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
