LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep:
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas
	lmr prep-grid --job=grid.lua
	lmr prep-sim --job=transient.lua

run:
	lmr run

post:
	lmr snapshot2vtk --all --add-vars="mach,pitot"

test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
