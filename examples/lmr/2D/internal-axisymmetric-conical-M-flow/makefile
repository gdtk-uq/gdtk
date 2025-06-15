help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep: contour gas init

contour:
	python m_flow.py

gas:
	echo "model = 'ideal gas'; species = {'air'}" > ideal-air.lua
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas
	rm ideal-air.lua

init:
	lmr prep-grid
	lmr prep-sim

run:
	lmr run

post:
	lmr snapshot2vtk --add-vars="mach"

wall-data:
	lmr slice-flow --add-vars="mach" --slice-list="0,:,49,0" --output="wall-0.dat"
	lmr slice-flow --add-vars="mach" --slice-list="0,:,48,0" --output="wall-1.dat"
	lmr slice-flow --add-vars="mach" --slice-list="0,:,47,0" --output="wall-2.dat"
	lmr slice-flow --add-vars="mach" --slice-list="0,:,46,0" --output="wall-3.dat"

