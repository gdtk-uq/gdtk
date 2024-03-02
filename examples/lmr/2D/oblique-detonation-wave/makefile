LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep: job.lua
	lmr prep-grid
	lmr prep-sim

run:
	lmr run

post: vtk

vtk:
	lmr snapshot2vtk --all

clean:
	- rm -rf $(LMR_DIR)
