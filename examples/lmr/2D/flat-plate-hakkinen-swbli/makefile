LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run-transient"
	@echo -e "   post"

prep:
	prep-gas ideal-air.lua ideal-air.gas
	lmr prep-grid --job=grid.lua

run-transient:
	lmr prep-sim --job=transient.lua
	lmr run

post:
	lmr snapshot2vtk --all
	lmr slice-flow --output=surface.data \
        --slice-list="0,:,0,0;2,:,0,0;4,:,0,0;6,:,0,0;8,:,0,0;10,:,0,0;12,:,0,0;14,:,0,0;16,:,0,0;18,:,0,0"

clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
