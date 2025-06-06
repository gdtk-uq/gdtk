LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep:
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/air-5sp-1T.inp ./air-5sp.inp
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua .
	lmr prep-gas -i air-5sp.inp -o air-5sp.lua
	lmr prep-chem -g air-5sp.lua -i GuptaEtAl-air-reactions.lua -o air-5sp-6r.lua
	lmr prep-grid
	lmr prep-sim

run:
	mpirun -np 4 lmr-mpi-run | tee log.txt

post:
	lmr snapshot2vtk --all

standoff:
	lmr custom-script --job=shock-standoff.lua

clean:
	- rm -rf $(LMR_DIR) air-5sp.inp air-5sp.lua GuptaEtAl-air-reactions.lua air-5sp-6r.lua log.txt
