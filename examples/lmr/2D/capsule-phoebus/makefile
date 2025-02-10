FILES := lmrsim ideal-air.gas

gas: ideal-air.inp
	lmr prep-gas -i ideal-air.inp -o ideal-air.gas

grid: grid.lua
	lmr prep-grids --job=grid.lua --verbose

flow: sim.lua
	lmr prep-flow --job=sim.lua --verbose

prep: gas grid flow

run-shared:
	lmr run

run:
	mpirun --use-hwthread-cpus -np 8 lmrZ-mpi-run

plot-diagnostics:
	lmr plot-diagnostics --style=kad --mode=live

post:
	lmr snapshot2vtk --all --add-vars="Mach,total-h"

plot-residuals:
	lmr residual2vtk --all

clean:
	rm -rf $(FILES)
