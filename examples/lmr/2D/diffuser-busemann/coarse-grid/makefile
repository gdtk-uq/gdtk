gas:
	printf '%s\n' 'model = "ideal gas"' 'species = {"air"}' > ideal-air.inp
	lmr prep-gas -i ideal-air.inp -o ideal-air.lua

grid:
	./generate-contour.py
	sed -i '1 s/./# &/' bd-contour.dat
	lmr prep-grid

prep-sim: gas
	lmr prep-sim

prep: gas grid prep-sim

run:
	mpirun -np 4 lmrZ-mpi-run

vtk:
	lmr snapshot2vtk --add-vars="mach,total-p"

wall-props:
	lmr slice-flow -l ":,:,$$,0" --output=wall-properties.dat

clean:
	- rm -rf lmrsim
	- rm ideal-air.lua
	- rm bd-contour.dat bd-props.txt
	- rm wall-properties.dat

