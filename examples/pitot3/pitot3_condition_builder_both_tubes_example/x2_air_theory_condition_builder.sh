# x2_air_theory_condition_builder.sh: 
#
# This is a shell script to run the PITOT3 condition builder.
# This is an example with a single driver condition, a linearly spaced shock tube fill pressure, and a logarithmically spaced acceleration tube fill pressure.
# Obviously either tube can be linear or log, but this is just an example showing they can be different.
# Chris James (c.james4@uq.edu.au) - 15/03/24

pitot3_condition_builder.py --config_file x2_air_theory_condition_builder.yaml
