
nope:
	echo "nope"

clean:
	rm -rf lmrsim ideal-air.lua
			
post:
	lmr snapshot2vtk -a

prep:
	lmr prep-gas -i ideal-air.inp -o ideal-air.lua
	lmr prep-grid --job=grid.lua
	lmr prep-flow --job=elec.lua

run:
	lmr-run

mes:
	python3 mes.py
