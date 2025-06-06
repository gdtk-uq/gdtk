LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep:
	lmr prep-gas -i tp-air.inp -o tp-air.lua
	lmr prep-grid
	lmr prep-sim

run:
	mpirun -np 4 lmr-mpi-run | tee log.txt

post:
	lmr snapshot2vtk --all

standoff:
	lmr custom-script --job=shock-standoff.lua

clean:
	- rm -rf $(LMR_DIR) tp-air.lua log.txt
