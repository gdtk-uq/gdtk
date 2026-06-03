
nope:
	echo "Nope"

clean:
	rm -rf tp-air.lua lmrsim
			
prep:
	lmr prep-gas -i tp-air.inp -o tp-air.lua
	lmr prep-grid
	lmr prep-sim

post:
	lmr snapshot2vtk -a

run:
	mpirun -np 14 lmrZ-mpi-run

