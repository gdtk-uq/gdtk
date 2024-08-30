LMR_DIR := lmrsim 

help:
	echo "Look in the makefile..."

grid:
	ugrid_partition sphere.su2 mapped_cells 1 2
	mkdir su2_grid
	mv block_* su2_grid/
	mv mapped_cells su2_grid/

prep:
	prep-gas gm-air-7sp-2T.inp air-7sp-2T.gas
	prep-chem air-7sp-2T.gas rm-park-7sp-2T.inp air-7sp-2T.chem
	prep-kinetics air-7sp-2T.gas air-7sp-2T.chem eem-park-7sp-2T.inp air-7sp-2T.exch
	lmr prep-grid --with-cea-usage
	lmr prep-sim --with-cea-usage --job=job.lua

run:
	/usr/bin/time -f "RSS=%M elapsed=%E" mpirun -np 1 lmrZ-mpi-run | tee log.dat

check:
	lmrZ-check-jacobian --read-frozen-limiter-values=true --output=norms.dat


post: vtk

vtk:
	lmr snapshot2vtk -s 5 --add-vars=mach

clean:
	- rm -rf $(LMR_DIR)
	- rm *.gas *.chem *.exch
	- rm thermo.inp thermo.lib thermo.out trans.inp trans.lib trans.out tmp.inp tmp.out
	- rm log.dat
	- rm *png
	- rm -r su2_grid
