
nope:
	echo "Nope"

grid:
	lmr custom-script --job=gengrid.lua

clean:
	rm -rf lmrsim ideal-air.lua

deepclean: clean
	rm -rf grid.su2 su2grid mapped_cells

partition:
	ugrid_partition grid.su2 mapped_cells 8 2
	mkdir -p su2grid
	mv block_* su2grid/

run:
	mpirun -np 8 lmrZ-mpi-run
			
prep:
	prep-gas ideal-air.inp ideal-air.lua
	lmr prep-grid --job=fp-grid.lua
	lmr prep-flow --job=fp.lua

post:
	lmr snapshot2vtk

