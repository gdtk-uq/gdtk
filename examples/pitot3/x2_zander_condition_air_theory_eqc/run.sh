# x2_zander_condition_air_theory_custom_driver_condition.sh: 
# Shell script to run the Zander condiiton
# air theoretical example for PITOT3 with a custom driver condition specified.
# Chris James (c.james4@uq.edu.au) - 31/08/21

prep-gas tp-he-ar.inp  tp-he-ar.lua
prep-gas tp-n2-o2-ions.inp tp-n2-o2-ions.lua
pitot3.py --config_file input.yaml
