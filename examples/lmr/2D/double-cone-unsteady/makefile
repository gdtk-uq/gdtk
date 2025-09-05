LMR_DIR := lmrsim 

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep: gas grid init

gas:
	prep-gas gm-ideal-N2.inp ideal-N2.gas

grid:
	lmr custom-script --job=grid.lua
	ugrid_partition grid.su2 mapped_cells 8 2 true
	mkdir su2grid
	mv grid.su2 block* mapped_cells su2grid
	lmr prep-grid

init:
	lmr prep-sim --job=job.lua

run:
	mpirun -np 8 lmrZ-mpi-run -v | tee log.dat

post: vtk

vtk:
	lmr snapshot2vtk --all --add-vars=mach

clean:
	- rm -rf $(LMR_DIR)
	- rm *.gas *.chem *.exch
	- rm log.dat
	- rm *png
	- rm -r su2grid
