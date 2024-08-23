
gas:
	lmr prep-gas -i gm-air11-2T.lua -o air11-2T.gas

chem:
	lmr prep-chem -g air11-2T.gas -i rr-park-air11-2T.lua -o park-air11-2T.chem

kinetics:
	lmr prep-energy-exchange -g air11-2T.gas -i ee-park-air11-2T.lua -o park-air11-2T.exch

grid:
	lmr prep-grid -j grid.lua

prep-sim:
	lmr prep-sim

prep: gas chem kinetics grid prep-sim

run:
	mpirun -np 16 lmr-mpi-run | tee LOGFILE

post:
	lmr snapshot2vtk -f --add-vars="mach"
