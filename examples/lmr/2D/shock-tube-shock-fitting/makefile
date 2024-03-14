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
	python3 sod-driving-left.py > sod-driving-left.transcript
	lmr slice-flow --output=eilmer.data --slice-list=":,:,0,0"
	./plot.sh

clean:
	- rm -rf $(LMR_DIR)
	- rm *.eps eilmer.data analytic.data sod-driving-left.transcript
	- rm ideal-air.gas
