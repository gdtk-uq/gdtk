LMR_DIR := lmrsim 

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep: partition gas grid init

gas:
	prep-gas gm-ideal-air.inp ideal-air.gas

grid:
	lmr prep-grid

init:
	lmr prep-sim --job=job.lua

partition:
	ugrid_partition ramp.su2 mapped_cells 1 2 true

run:
	/usr/bin/time -f "RSS=%M elapsed=%E" mpirun -np 1 lmrZ-mpi-run | tee log.dat

post: vtk

vtk:
	lmr snapshot2vtk --all --add-vars=mach
	lmr residual2vtk --all
	lmr limiter2vtk --all
	lmr gradient2vtk --all

clean:
	- rm -rf $(LMR_DIR)
	- rm *.gas *.chem *.exch
	- rm log.dat
	- rm block* mapped_cells
