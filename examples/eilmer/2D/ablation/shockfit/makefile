filename=sc

nope:
	echo "Nope"

shockfit:
	python3 inflate_sf_shock.py 4e-4

clean:
	rm -rf flow grid hist loads plot solid config *.log uniform.o* uniform.e* gm-air11-2T.lua rr-kim-air11-2T.lua ee-kim-air11-2T.lua log*.txt
			
prep:
	prep-gas gm-air11-2T.inp gm-air11-2T.lua
	prep-chem gm-air11-2T.lua rr-kim-air11-2T.inp rr-kim-air11-2T.lua
	prep-kinetics gm-air11-2T.lua rr-kim-air11-2T.lua ee-kim-air11-2T.inp ee-kim-air11-2T.lua
	e4shared --prep --job=${filename}

run:
	mpirun -np 8 e4mpi --run --job=${filename} --tindx-start=last --report-residuals | tee -a log.txt

post:
	e4shared --post --job=${filename} --vtk-xml --tindx-plot=all
