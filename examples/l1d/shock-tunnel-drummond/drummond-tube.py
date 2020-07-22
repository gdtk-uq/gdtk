# drummond-tube.py
#
# Tube wall definition for the Drummond Shock Tunnel.
#
tube.n = 4000
add_break_point(-3.785,  0.0585)    # upstream-end of the driver tube
add_break_point(-3.035,  0.0585)
add_break_point(-3.015,  0.0620)    # steel-diaphragm station
add_break_point( 0.000,  0.0620)    # downstrem end of shock tube
add_break_point( 0.043,  0.0620)    # start of contraction to throat
add_break_point( 0.080,  0.0220)    # start of throat
add_break_point( 0.100,  0.0220)    # start of nozzle (conical shape)
add_break_point( 0.2653, 0.0700)    # start of parallel test section
add_break_point( 0.30,   0.0700)
#
add_history_loc(-0.295) # 0, heat flux gauge
add_history_loc(-0.078) # 1, pressure transducer
add_history_loc( 0.000) # 2, joint at nozzle block
add_history_loc( 0.090) # 3, mid-point of nozzle throat
add_history_loc( 0.265) # 4, nozzle exit plane
