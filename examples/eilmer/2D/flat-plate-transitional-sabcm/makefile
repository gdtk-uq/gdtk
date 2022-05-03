filename=fp

all: clean prep run post

clean:
	rm -rf residuals flow grid hist loads plot solid config log.txt ideal-air-gas-model.lua ${filename}-ref-residuals.saved e4-nk.diagnostics.dat wht_*.svg

prep:
	prep-gas ideal-air.inp ideal-air-gas-model.lua
	e4shared --prep --job=${filename}

run:
	mpirun -np 8 e4-nk-dist --job=${filename} | tee log.txt

post:
	e4shared --post --job=${filename} --vtk-xml --tindx-plot=all
