LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep: grid init-sim

grid: grid.lua
	lmr prep-grid --job=grid.lua

init-sim: let.lua
	lmr prep-sim --job=let.lua

run:
	# lmr run
	mpirun -np 7 --oversubscribe lmr-mpi-run

post: vtk slice

vtk:
	lmr snapshot2vtk --all

slice:
	lmr slice-flow --slice-list=":,:,0,0"

clean:
	- rm -rf $(LMR_DIR)
