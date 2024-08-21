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

post: vtk probe slice

vtk:
	lmr snapshot2vtk --all

probe:
	lmr probe-flow --location=0.28,0,0 --names=pos.x,T,vel.x,p

slice:
	lmr slice-flow --slice-list=0,:,0,0 --output=profile.dat --names=T,vel.x,p

plots:
	gnuplot plot-profile.gplot

test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
	- rm profile.dat *-profile.pdf
