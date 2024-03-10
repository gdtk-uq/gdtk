
nope:
	echo "Nope"

clean:
	rm -rf lmrsim ideal-air.gas

run:
	mpirun -np 8 lmrZ-mpi-run
			
prep:
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas
	lmr prep-grid
	lmr prep-flow

post:
	lmr snapshot2vtk

