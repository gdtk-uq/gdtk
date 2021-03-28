filename=fire

clean:
	rm -rf flow grid hist loads plot solid config *.log fire_*.o* fire_*.e*
			
prep:
	prep-gas air-11sp-gas-model.inp air-11sp-gas-model.lua
	prep-chem air-11sp-gas-model.lua air-11sp-reactions.inp air-11sp-reactions.lua
	e4shared --prep --job=${filename}

run:
	mpirun -np 16 e4mpi --run --job=${filename} --tindx-start=0

post:
	e4shared --post --job=${filename} --vtk-xml --tindx-plot=all
