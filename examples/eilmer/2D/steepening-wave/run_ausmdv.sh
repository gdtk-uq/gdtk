#!/bin/bash
# run.sh

rm -rf ausmdv
mkdir -p ausmdv

cp -r blank ausmdv/032
sed -i 's/N=128/N=32/' ausmdv/032/swp.lua
sed -i 's/asf/ausmdv/' ausmdv/032/swp.lua

cp -r blank ausmdv/064
sed -i 's/N=128/N=64/' ausmdv/064/swp.lua
sed -i 's/asf/ausmdv/' ausmdv/064/swp.lua

cp -r blank ausmdv/128
sed -i 's/N=128/N=128/' ausmdv/128/swp.lua
sed -i 's/asf/ausmdv/' ausmdv/128/swp.lua


cd ausmdv/032
prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=swp
e4shared --run --job=swp --verbosity=1
cd ../..

cd ausmdv/064
prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=swp
e4shared --run --job=swp --verbosity=1
cd ../..

cd ausmdv/128
prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=swp
e4shared --run --job=swp --verbosity=1
cd ../..

python3 scripts/post_process.py ausmdv/032 ausmdv/064 ausmdv/128
