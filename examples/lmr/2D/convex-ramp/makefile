LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep:
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-5sp-gas-model.lua .
	lmr prep-gas -i air-5sp-gas-model.lua -o air-5sp-2T.gas
	lmr prep-grid --job=grid.lua

prep-transient: transient.lua
	lmr prep-sim --job=transient.lua

prep-steady: steady.lua
	lmr prep-sim --job=steady.lua

run-transient:
	lmr run

run-steady:
	mpirun --use-hwthread-cpus -np 8 lmrZ-mpi-run

post:
	lmr snapshot2vtk --all

test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v

clean:
	- rm -rf $(LMR_DIR)
	- rm air-5sp-gas-model.lua air-5sp-2T.gas
