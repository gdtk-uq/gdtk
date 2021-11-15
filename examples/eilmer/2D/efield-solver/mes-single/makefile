
nope:
	echo "nope"

clean:
	rm -rf flow grid hist loads plot solid config CellData ideal-air-gas-model.lua npdata
			
post:
	e4shared --post --job=elec --vtk-xml --tindx-plot=1 --plotTag="efield"

prep:
	prep-gas ideal-air.inp ideal-air-gas-model.lua
	e4shared --prep --job=elec

run:
	e4shared --run --job=elec --max-cpus=1

mes:
	python3 mes.py
