
nope:
	echo "Nope"

clean:
	rm -rf flow grid hist loads plot solid config limiter-values residuals swbli-ref-residuals.saved ideal-air-gas-model.lua log.txt gm-*.lua rr-*.lua e4-nk.diagnostics.dat

run:
	mpirun -np 16 e4-nk-dist --job=swbli | tee -a log.txt
			
			
prep:
	prep-gas gm-air5.inp gm-air5.lua
	prep-chem gm-air5.lua rr-kim-air5.inp rr-kim-air5.lua
	e4shared --prep --job=swbli

post:
	e4shared --post --job=swbli --vtk-xml --tindx-plot=last

