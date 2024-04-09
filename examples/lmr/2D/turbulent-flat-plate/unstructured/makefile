
nope:
	echo "Nope"

grid:
	lmr custom-script --job=gengrid.lua

clean:
	rm -rf lmrsim

deepclean: clean
	rm -rf grid.su2 su2grid mapped-cells

partition:
	ugrid_partition grid.su2 mapped-cells 8 2
	mkdir -p su2grid
	mv block_* su2grid/

run:
	mpirun -np 8 lmrZ-mpi-run
			
prep:
	mkdir -p lmrsim
	cp mapped-cells lmrsim
	lmr prep-gas -i ideal-air.lua -o lmrsim/ideal-air.gas
	lmr prep-grid --job=fp-grid.lua
	lmr prep-flow --job=fp.lua

post:
	lmr snapshot2vtk

