
nope:
	echo "Nope"

grid:
	lmr custom-script --job=gengrid.lua

clean:
	rm -rf lmrsim

deepclean: clean
	rm -rf mapped-cells grid.su2 su2grid

partition:
	ugrid_partition grid.su2 mapped-cells 4 2
	mkdir -p su2grid
	mv block_* su2grid/

prep: job.lua
	mkdir -p lmrsim
	cp mapped-cells lmrsim/mapped-cells
	lmr prep-gas -i species.inp -o lmrsim/species.lua
	lmr prep-chem -g lmrsim/species.lua -i reactions.inp -o lmrsim/reactions.lua
	lmr prep-grid -j importgrid.lua
	lmr prep-sim

run:
	mpirun -np 4 lmrZ-mpi-run

transient:
	mpirun -np 4 lmr-mpi-run

post:
	lmr snapshot2vtk --all
