filename=sc

nope:
	echo "Nope"

shockfit:
	python3 ../../../scripts/inflate_sf_shock.py ../parameters.py

clean:
	rm -rf flow grid hist loads plot solid config *.log gm-air11-2T.lua rr-kim-air11-2T.lua ee-kim-air11-2T.lua log*.txt *shock_shape.dat sc-ref-residuals.saved e4-nk.diagnostics.dat residuals
			
prep:
	prep-gas gm-air11-2T.inp gm-air11-2T.lua
	prep-chem gm-air11-2T.lua rr-kim-air11-2T.inp rr-kim-air11-2T.lua
	prep-kinetics gm-air11-2T.lua rr-kim-air11-2T.lua ee-kim-air11-2T.inp ee-kim-air11-2T.lua
	e4shared --prep --job=${filename}

run:
	mpirun -np 8 e4mpi --job=${filename} --run --tindx-start=last | tee -a log.txt

post:
	e4shared --post --job=${filename} --vtk-xml --tindx-plot=all
