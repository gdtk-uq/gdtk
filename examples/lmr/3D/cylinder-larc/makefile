LMR_DIR := lmrsim 

help:
	echo "Look in the makefile..."

grid:
	ugrid_partition cylinder.su2 mapped_cells 8 3 true
	mkdir su2_grid
	mv block_* su2_grid/
	mv mapped_cells su2_grid/

prep:
	prep-gas gm-air-5sp-2T.inp air-5sp-2T.gas
	prep-chem air-5sp-2T.gas rm-park-5sp-2T.inp air-5sp-2T.chem
	prep-kinetics air-5sp-2T.gas air-5sp-2T.chem eem-park-5sp-2T.inp air-5sp-2T.exch
	lmr prep-grid --with-cea-usage
	lmr prep-sim --with-cea-usage --job=job.lua

run:
	/usr/bin/time -f "RSS=%M elapsed=%E" mpirun -np 8 lmrZ-mpi-run | tee log.dat

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
