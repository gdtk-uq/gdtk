# t4-tube.py
#
# Tube wall definition for T4 with Mach 7 nozzle, as at 2019-03-21.
# From L1d3 simulation of shot 11315
#
# Note that, where ever a restriction is located, we need to have a significant
# length of tube for the gas cells to actually see the correct tube area.
# Also, the transitions between constant-area sections need to be smooth
# enough for the gas-dynamic simulation to be well behaved.
tube.n = 10000 # number of linear segments used for fast interpolation
add_break_point(-4.87, 0.442)  # upstream-end of new compressed-air reservoir
add_break_point(-0.90, 0.442)
add_break_point(-0.80, 0.408)  # start of sleeve
add_break_point(-0.30, 0.408)  # end of sleeve
add_break_point(-0.20, 0.238)  # start of equivalent duct modelling slots
add_break_point(-0.10, 0.238)  # end of slots.
add_break_point( 0.00, 0.229)  # upstream-end of compression tube
add_break_point(25.91, 0.229)  # downstream-end of compression tube
add_break_point(26.00, 0.076)  # full bore of shock tube
add_break_point(36.06, 0.076)  # downstream-end of shock tube (36.045 + T4 M7 nozzle 0.012)
add_break_point(36.08, 0.021)  # start of nozzle throat (T4 M7 nozzle)
add_break_point(36.10, 0.021)  # end of nozzle throat (T4 M7 nozzle)
add_break_point(37.00, 0.262)  # exit-plane of nozzle, start of test section
#
# Sensor locations...
# Remember that these history locations will be numbered from 0 in the set
# of history data files.  Thus, the nozzle-supply sensor will be numbered 4.
add_history_loc(25.8)   # downstream-end of the compression tube
add_history_loc(30.01)   # shock-tube station 1
add_history_loc(32.00)   # shock-tube station 2
add_history_loc(34.01)   # shock-tube station 3
add_history_loc(36.01)   # downstream-end of shock tube (nozzle-supply region)
add_history_loc(37.00)   # downstream-end of nozzle
