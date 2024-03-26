LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

prep:
	cp ${DGD_REPO}/examples/kinetics/hydrogen-ignition-delay/Rogers-Schexnayder-species.inp ./species.lua
	cp ${DGD_REPO}/examples/kinetics/hydrogen-ignition-delay/Rogers-Schexnayder.lua ./reactions.lua
	sed -i '/^options.H2_O2_only =/s/false/true/' species.lua
	sed -i '/^options.H2_O2_only =/s/false/true/' reactions.lua
	lmr prep-gas -i species.lua -o h2o2.gas
	lmr prep-reactions -g h2o2.gas -i reactions.lua -o h2o2.chem
	lmr prep-grid --job=grid.lua
	lmr prep-sim --job=transient.lua

run:
	lmr run

post:
	lmr snapshot2vtk --all

test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v

clean:
	- rm -rf $(LMR_DIR)
	- rm species.lua h2o2.gas reactions.lua h2o2.chem
