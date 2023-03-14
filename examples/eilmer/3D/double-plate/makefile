
nope:
	echo "Nope"

grid: gengrid.lua
	e4shared --custom-script --script-file="gengrid.lua"

deepclean: clean
	rm -rf grid.su2 su2grid mapped_cells

clean:
	rm -rf flow grid hist loads plot solid config limiter-values residuals dp-ref-residuals.saved ideal-air.lua log.txt gm-*.lua rr-*.lua e4-nk.diagnostics.dat CellData

run:
	mpirun -np 16 e4-nk-dist --job=dp --snapshot-start=last | tee -a log.txt
			
prep:
	prep-gas ideal-air.inp ideal-air.lua
	e4shared --prep --job=dp

post:
	e4shared --post --job=dp --vtk-xml --tindx-plot=last

residual:
	e4shared --post --job=dp --vtk-xml --tindx-plot=last --plotTag="residual"

partition:
	ugrid_partition grid.su2 mapped_cells 16 3
	mkdir -p su2grid
	mv block_* su2grid/

