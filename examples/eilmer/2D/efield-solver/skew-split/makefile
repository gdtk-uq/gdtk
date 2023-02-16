
nope:
	echo "nope"

clean:
	rm -rf flow grid hist loads plot solid config CellData ideal-air-gas-model.lua npdata ideal-air.lua
			
post:
	e4shared --post --job=elec --vtk-xml --tindx-plot=1 --plotTag="efield"

prep:
	prep-gas ideal-air.inp ideal-air.lua
	e4shared --prep --job=elec

run:
	e4shared --run --job=elec --max-cpus=1

mpirun:
	mpirun -np 2 e4mpi --run --job=elec

mes:
	mkdir -p npdata
	apply_filter.sh ../scripts/mes.py plot/elec.efield-t0001.pvtu
