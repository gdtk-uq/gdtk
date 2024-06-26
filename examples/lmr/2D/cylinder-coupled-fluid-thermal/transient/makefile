LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "	prep"
	@echo "	run/mpirun"
	@echo "	post"

prep: prep-gas grid init

prep-gas: ideal-air.lua
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas

grid:
	lmr prep-grid

init:
	lmr prep-sim

run:
	lmr run

mpirun:
	mpirun -np 4 lmr-mpi-run

post: vtk

vtk:
	lmr snapshot2vtk --all

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
