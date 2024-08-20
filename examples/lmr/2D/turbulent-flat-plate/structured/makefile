
nope:
	echo "Nope"

clean:
	rm -rf lmrsim ideal-air.gas test.lua __pycache__

run:
	mpirun -np 8 lmrZ-mpi-run
			
prep:
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas
	lmr prep-grid --job=job.lua
	lmr prep-flow --job=job.lua

post:
	lmr snapshot2vtk

