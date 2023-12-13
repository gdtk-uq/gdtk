#!/bin/bash
# run.sh

rm -rf ldfss2
mkdir -p ldfss2

cp -r blank ldfss2/032
sed -i 's/N=128/N=32/' ldfss2/032/swp.lua
sed -i 's/asf/ldfss2/' ldfss2/032/swp.lua

cp -r blank ldfss2/064
sed -i 's/N=128/N=64/' ldfss2/064/swp.lua
sed -i 's/asf/ldfss2/' ldfss2/064/swp.lua

cp -r blank ldfss2/128
sed -i 's/N=128/N=128/' ldfss2/128/swp.lua
sed -i 's/asf/ldfss2/' ldfss2/128/swp.lua


cd ldfss2/032
prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=swp
e4shared --run --job=swp --verbosity=1
cd ../..

cd ldfss2/064
prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=swp
e4shared --run --job=swp --verbosity=1
cd ../..

cd ldfss2/128
prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=swp
e4shared --run --job=swp --verbosity=1
cd ../..

python3 scripts/post_process.py ldfss2/032 ldfss2/064 ldfss2/128
