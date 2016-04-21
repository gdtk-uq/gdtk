python make_source_terms.py
e4shared --job=mms --prep
e4shared --job=mms --run
e4shared --job=mms --post --tindx-plot=20 --ref-soln=udf-bc.lua  --norms="rho" | tail -3 > rho-norms.txt
