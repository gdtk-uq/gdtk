LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep-A"
	@echo "   run"
	@echo "   post"

prep-A: ideal-air.lua grid.lua transient-A.lua
	prep-gas ideal-air.lua ideal-air.gas
	python3 prepare-static-profile.py > prepare-static-profile.transcript
	lmr prep-grid --job=grid.lua
	lmr prep-sim --job=transient-A.lua

run:
	lmr run

post:
	lmr snapshot2vtk --all
	python3 classic-shock-tube.py > classic-shock-tube.transcript
	lmr slice-flow --final --output=eilmer.data --slice-list=:,:,0,0
	./plot.sh

clean:
	- rm -rf $(LMR_DIR)
	- rm static-profile.data prepare-static-profile.transcript
	- rm transient-profile.zip prepare-transient-profile.transcript
	- rm *.eps eilmer.data analytic.data classic-shock-tube.transcript
	- rm ideal-air.gas
