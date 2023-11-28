LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"

prep: prep-gas grid

prep-gas: gm-ideal-air.inp
	prep-gas gm-ideal-air.inp ideal-air.gas

grid:
	lmr prep-grid

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas

