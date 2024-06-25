
nope:
	echo "Nope"

grid: gengrid.lua
	lmr custom-script --job="gengrid.lua"

partition:
	ugrid_partition grid.su2 mapped-cells 16 2
	mkdir -p su2grid
	mv block_* su2grid/

deepclean: clean
	rm -rf mapped-cells grid.su2 su2grid

clean:
	rm -rf lmrsim

run:
	mpirun -np 16 lmrZ-mpi-run
			
prep:
	mkdir -p lmrsim
	cp mapped-cells lmrsim
	lmr prep-gas -i gm-air.inp -o lmrsim/gm-air.lua
	lmr prep-grid --job=af-grid.lua
	lmr prep-flow --job=af.lua

post:
	lmr snapshot2vtk

postall:
	lmr snapshot2vtk --all

all: deepclean grid partition prep run post
