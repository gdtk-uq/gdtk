LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   grid"

prep: prep-gas grid

prep-gas: ideal-air.lua
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas

grid:
	lmr prep-grid

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas

