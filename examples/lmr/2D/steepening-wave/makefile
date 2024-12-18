
nope:
	echo "Nope"

clean:
	rm -rf asf ausmdv

prep:
	mkdir -p asf
	cp -r blank asf/032
	sed -i 's/N=128/N=32/' asf/032/swp.lua
	sed -i 's/N=128/N=32/' asf/032/swp-grid.lua
	cd asf/032; lmr prep-gas -i ideal-air.lua -o ideal-air.gas; lmr prep-grid -j swp-grid.lua; lmr prep-sim -j swp.lua
	
	cp -r blank asf/064
	sed -i 's/N=128/N=64/' asf/064/swp.lua
	sed -i 's/N=128/N=64/' asf/064/swp-grid.lua
	cd asf/064; lmr prep-gas -i ideal-air.lua -o ideal-air.gas; lmr prep-grid -j swp-grid.lua; lmr prep-sim -j swp.lua
	
	cp -r blank asf/128
	sed -i 's/N=128/N=128/' asf/128/swp.lua
	sed -i 's/N=128/N=128/' asf/128/swp-grid.lua
	cd asf/128; lmr prep-gas -i ideal-air.lua -o ideal-air.gas; lmr prep-grid -j swp-grid.lua; lmr prep-sim -j swp.lua
			
run:
	cd asf/032; lmr run
	cd asf/064; lmr run
	cd asf/128; lmr run

post:
	python3 scripts/post_process.py asf/032 asf/064 asf/128
