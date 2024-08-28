LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep: gas grid init-transient

gas: ideal-air.lua
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas

grid: grid.lua
	lmr prep-grid --job=grid.lua

init-transient: transient.lua
	lmr prep-sim --job=transient.lua

run:
	lmr run

post: vtk trajectory-plot

vtk:
	lmr snapshot2vtk --all

trajectory-plot:
	ruby analytic-trajectory.rb > analytic-trajectory.data
	gnuplot plot-trajectory.gnuplot
test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas userPadSave.lua
	- rm piston.data analytic-trajectory.data
	- rm piston-position.pdf piston-velocity.pdf
