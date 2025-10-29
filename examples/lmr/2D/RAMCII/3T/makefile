prep:
	lmr prep-gas -i gas.lua -o gas.gas | tee log.txt
	lmr prep-chem -g gas.lua -i chem.lua -o chem.chem | tee -a log.txt
	lmr prep-kinetics -g gas.gas -r chem.chem -i kinetics.lua -o kinetics.kin | tee -a log.txt
	lmr prep-grid --job=grid.lua | tee -a log.txt
	lmr prep-sim --job=ram.lua | tee -a log.txt

clean:
	-rm -f gas.gas kinetics.kin chem.chem
	-rm -r lmrsim
	-rm log.txt

run:
	mpirun -n 8 lmrZ-mpi-run | tee -a log.txt

post:
	lmr snapshot2vtk -a --add-vars="mach,total-h,molef,conc"
