filename=cat

all: clean prep run post

clean:
	rm -rf flow grid hist loads plot solid config *.log catalytic.o* catalytic.e* air-5sp-gas-model.lua air-5sp-reactions.lua
			
prep:
	prep-gas air-5sp-gas-model.inp air-5sp-gas-model.lua
	prep-chem air-5sp-gas-model.lua air-5sp-reactions.inp air-5sp-reactions.lua
	e4shared --prep --job=${filename}

run:
	mpirun -np 16 e4mpi --run --job=${filename}

post:
	e4shared --post --job=${filename} --vtk-xml --tindx-plot=all
