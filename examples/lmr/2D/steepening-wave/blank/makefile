
nope:
	echo "Nope"

clean:
	rm -rf lmrsim ideal-air.gas

run:
	lmr run
			
prep:
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas
	lmr prep-grid -j swp-grid.lua
	lmr prep-sim -j swp.lua

post:
	lmr snapshot2vtk

