python make_lua_files.py
e4shared --job=mms --prep
e4shared --job=mms --run
e4shared --job=mms --post --tindx-plot=last --ref-soln=ref-soln.lua  --norms="rho" | tail -3 > rho-norms.txt
