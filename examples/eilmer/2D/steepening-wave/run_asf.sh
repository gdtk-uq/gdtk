#!/bin/bash
# run.sh

rm -rf asf
mkdir -p asf

cp -r blank asf/032
sed -i 's/N=128/N=32/' asf/032/swp.lua

cp -r blank asf/064
sed -i 's/N=128/N=64/' asf/064/swp.lua

cp -r blank asf/128
sed -i 's/N=128/N=128/' asf/128/swp.lua


cd asf/032
prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=swp
e4shared --run --job=swp --verbosity=1
cd ../..

cd asf/064
prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=swp
e4shared --run --job=swp --verbosity=1
cd ../..

cd asf/128
prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=swp
e4shared --run --job=swp --verbosity=1
cd ../..

python3 scripts/post_process.py asf/032 asf/064 asf/128
