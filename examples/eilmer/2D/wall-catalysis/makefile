filename=cat

all: clean prep run post

clean:
	rm -rf flow grid hist loads plot solid config *.log catalytic.o* catalytic.e* \
		air-5sp-gas-model.lua GuptaEtAl-air-reactions-2T.lua air-energy-exchange.lua gm.lua rr.lua ee.lua

prep:
	./prep-gas-and-chem.sh
	e4shared --prep --job=${filename}

run:
	mpirun -np 16 --oversubscribe e4mpi --run --job=${filename}

post:
	e4shared --post --job=${filename} --vtk-xml --tindx-plot=all
