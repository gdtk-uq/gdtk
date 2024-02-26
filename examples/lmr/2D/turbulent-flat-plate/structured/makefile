
nope:
	echo "Nope"

clean:
	rm -rf lmrsim ideal-air.lua

run:
	mpirun -np 8 lmrZ-mpi-run
			
prep:
	prep-gas ideal-air.inp ideal-air.lua
	lmr prep-grid
	lmr prep-flow

post:
	lmr snapshot2vtk

