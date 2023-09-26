
nope:
	echo "Nope"

clean:
	rm -rf flow grid hist loads plot solid config limiter-values residuals swbli-ref-residuals.saved gm-*.lua rr-*.lua log.txt e4-nk.diagnostics.dat CellData

run:
	mpirun -np 16 e4-nk-dist --job=swbli | tee -a log.txt
			
			
prep:
	prep-gas gm-air.inp gm-air.lua
	e4shared --prep --job=swbli

post:
	e4shared --post --job=swbli --vtk-xml --tindx-plot=all

residual:
	e4shared --post --job=swbli --vtk-xml --tindx-plot=last --plotTag="residual"
