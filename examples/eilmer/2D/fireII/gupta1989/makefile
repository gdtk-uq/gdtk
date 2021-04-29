filename=fire

clean:
	rm -rf flow grid hist loads plot solid config log.txt fire_*.o* fire_*.e* air-11sp-*.lua
			
prep:
	prep-gas air-11sp-gas-model.inp air-11sp-gas-model.lua
	prep-chem air-11sp-gas-model.lua air-11sp-reactions.inp air-11sp-reactions.lua
	prep-kinetics air-11sp-gas-model.lua air-11sp-energy-exchange.inp air-11sp-energy-exchange.lua
	e4shared --prep --job=${filename}

run:
	mpirun -np 16 e4mpi --run --job=${filename} --tindx-start=0 | tee log.txt

post:
	e4shared --post --job=${filename} --vtk-xml --tindx-plot=all
