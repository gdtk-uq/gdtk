LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep: ideal-n2.lua grid.lua
	lmr prep-gas -i ideal-n2.lua -o ideal-n2.gas
	lmr prep-grid --job=grid.lua

run: transient.lua
	lmr prep-sim --job=transient.lua
	mpirun -np 7 --oversubscribe lmr-mpi-run

post:
	lmr snapshot2vtk --all --add-vars=mach,pitot

test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-n2.gas
