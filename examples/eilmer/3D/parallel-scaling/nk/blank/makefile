
np := $(shell lscpu --all --parse=CORE,SOCKET | grep -Ev "^\#" | sort -u | wc -l)
#OVERRIDE_CORE_COUNT

nope:
	echo "NOPE"

clean:
	rm -rf flow grid hist loads plot solid config residuals ideal-air-gas-model.lua e4-nk.diagnostics.dat bc-ref-residuals.saved log.txt

grid:
	e4shared --custom-script --script-file="gengrid.lua"

deepclean:
	rm -rf flow grid hist loads plot solid config residuals ideal-air-gas-model.lua e4-nk.diagnostics.dat bc-ref-residuals.saved log.txt su2grid mapped_cells grid.su2
			
prep: bc.lua
	prep-gas ideal-air.inp ideal-air-gas-model.lua
	e4shared --prep --job=bc

post:
	e4shared --post --job=bc --vtk-xml --tindx-plot=all

run:
	mpirun -np ${np} e4-nk-dist --job=bc --snapshot-start=last | tee -a log.txt

partition:
	ugrid_partition ../grid.su2 mapped_cells ${np} 3
	mkdir -p su2grid
	mv block_* su2grid/
