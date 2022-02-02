filename=bd

nope:
	echo "Nope"

clean:
	rm -rf flow grid hist loads plot solid config n2-*.lua log.txt profile.data
			
prep:
	prep-gas n2-3sp-gm.inp n2-3sp-gm.lua
	prep-chem n2-3sp-gm.lua n2-blank-rr.inp n2-blank-rr.lua
	prep-kinetics n2-3sp-gm.lua n2-blank-rr.lua n2-blank-ee.inp n2-blank-ee.lua
	e4shared --prep --job=${filename}

run:
	e4shared --job=${filename} --run | tee -a log.txt

post:
	e4shared --post --job=${filename} --vtk-xml --tindx-plot=all --add-vars="molef"

extract:
	e4shared --job=${filename} --post --slice-list="0,:,0,0" --output-file='profile.data'
