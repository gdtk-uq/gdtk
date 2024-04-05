LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep: grid.lua transient.lua
	lmr prep-gas -i two-temp-n2.lua -o two-temp-n2.gas
	lmr prep-chem -g two-temp-n2.gas -i chem.lua -o chem.chem
	lmr prep-energy-exchange -g two-temp-n2.gas -i VT-relaxation-time-selection.lua -o VT-relaxation-time-selection.kin
	lmr prep-grid --job=grid.lua
	lmr prep-sim --job=transient.lua


run:
	lmr run

post:
	lmr snapshot2vtk --all --add-vars=mach,pitot

clean:
	- rm -rf $(LMR_DIR)
	- rm two-temp-n2.gas 
	- rm VT-relaxation-time-selection.kin
	- rm chem.chem

