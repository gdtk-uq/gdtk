LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run-transient"
	@echo "   run-steady"
	@echo "   post"

prep:
	lmr prep-gas -i ideal-air.lua -o ideal-air.gas
	lmr prep-grid --job=grid.lua

run-transient:
	lmr prep-sim --job=transient.lua
	lmr run

run-steady:
	lmr prep-sim --job=steady.lua
	lmr run

post:
	lmr snapshot2vtk --all
	lmr slice-flow --output=surface.data \
        --slice-list="0,:,0,0;2,:,0,0;4,:,0,0;6,:,0,0;8,:,0,0;10,:,0,0;12,:,0,0;14,:,0,0;16,:,0,0;18,:,0,0"

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
