
nope:
	echo "Nope"

clean:
	rm -rf gm.lua rr.lua lmrsim
			
prep:
	lmr prep-gas -i gm-air6-1T.inp -o gm.lua
	lmr prep-chem -g gm.lua -i rr-jachimowski-air6.inp -o rr.lua
	lmr prep-grid
	lmr prep-sim

post:
	lmr snapshot2vtk -a

run:
	mpirun -np 16 lmrZ-mpi-run -s -1

nenzf1d:
	lmr prep-gas -i gm-air6-1T.inp -o gm.lua
	lmr prep-chem -g gm.lua -i rr-jachimowski-air6.inp -o rr.lua
	nenzf1d t4m7-air-11311.yaml
