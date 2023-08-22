python3 make_lua_files.py
e4shared --job=temporal-mms --prep
e4-nk-shared --job=temporal-mms --max-cpus=1
e4shared --job=temporal-mms --post --tindx-plot=last --ref-soln=ref-soln.lua  --norms="rho" --verbosity=0 | sed -n -e 6p -e 15p > rho-norms.txt
