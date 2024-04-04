LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep: grid.lua transient.lua
	lmr prep-grid --job=grid.lua
	lmr prep-sim --job=transient.lua


run:
	lmr run

post:
	lmr snapshot2vtk --all --add-vars=mach,pitot

clean:
	- rm -rf $(LMR_DIR)
