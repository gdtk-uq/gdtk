
copy-thermochem-files:
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-11sp-gas-model.lua .
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/Park-air-reactions-2T.lua .
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/Park-air-energy-exchange-2T.lua .
	
gas: copy-thermochem-files
	lmr prep-gas -i air-11sp-gas-model.lua -o air11-2T.gas

chem: gas
	lmr prep-chem -g air11-2T.gas -i Park-air-reactions-2T.lua -o park-air11-2T.chem

kinetics: gas
	lmr prep-energy-exchange -g air11-2T.gas -i Park-air-energy-exchange-2T.lua -o park-air11-2T.exch

grid:
	lmr prep-grid -j grid.lua

prep-sim:
	lmr prep-sim

prep: copy-thermochem-files gas chem kinetics grid prep-sim

run:
	mpirun -np 16 lmr-mpi-run | tee LOGFILE

post:
	lmr snapshot2vtk -f --add-vars="mach"

clean:
	- rm air-11sp-gas-model.lua
	- rm Park-air-reactions-2T.lua
	- rm Park-air-energy-exchange-2T.lua
	- rm -rf lmrsim/

