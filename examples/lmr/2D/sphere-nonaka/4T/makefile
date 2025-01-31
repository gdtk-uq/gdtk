LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep:
	lmr prep-gas -i air-5sp-4T-gas-model.lua -o air-5sp-4T.gas
	lmr prep-reactions -g air-5sp-4T.gas -i GuptaEtAl-air-reactions.lua -o air-5sp-6r-4T.chem
	lmr prep-energy-exchange -g air-5sp-4T.gas -r air-5sp-6r-4T.chem -i air-energy-exchange.lua -o air-energy-exchange-4T.exch
	lmr prep-grid
	lmr prep-sim

run:
	lmr run

post:
	lmr snapshot2vtk --all
	lmr custom-script --job=shock-shape.lua
	gnuplot plot-shock-shape.gplot
	python3 compute-error.py

slice:
	lmr slice-flow --final --slice-list="0,:,0,0" --names="T,T-N2,T-O2,T-NO" --output=stagnation-profile-T.dat
	lmr slice-flow --final --slice-list="0,:,0,0" --names="massf-N2,massf-O2,massf-N,massf-O,massf-NO" --output=stagnation-profile-species.dat

test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v

clean:
	- rm -rf $(LMR_DIR)
	- rm air-5sp-4T.gas air-5sp-6r-4T.chem air-energy-exchange-4T.exch
