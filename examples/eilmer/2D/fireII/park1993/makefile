filename=fire

clean:
	rm -rf flow grid hist loads plot solid config log.txt fire_*.o* fire_*.e* *-air11-2T.lua
			
prep:
	prep-gas gm-air11-2T.inp gm-air11-2T.lua
	prep-chem gm-air11-2T.lua rr-park-air11-2T.inp rr-park-air11-2T.lua
	prep-kinetics gm-air11-2T.lua ee-park-air11-2T.inp ee-park-air11-2T.lua
	e4shared --prep --job=${filename}

run:
	mpirun -np 16 e4mpi --run --job=${filename} --tindx-start=0 | tee log.txt

post:
	e4shared --post --job=${filename} --vtk-xml --tindx-plot=all
