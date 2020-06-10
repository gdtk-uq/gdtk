# x3r-tube.py
# X3R reflected shock tunnel facility as a 1-D tube model.
#
# Some history:
#   2020-06-02: Extracted by PJ from Sam Stennett's file for L1d3.
# More history:
#   Revised by David Gildfind, November 2012.
#   Edited by Andreas Andrianatos
#   Modified to X3R facility by Samuel Stennett
#     to account for the new reservoir extension

import math
import numpy as np

tube.n = 10000

# Component Lengths:
L_lead_in = 0.0500 # length of lead in on change in diameter
L_piston = 0.5000 # Length of piston
L_buffer = 0.1000 # Length of buffer
L_pp = 0.1000 # Length of pressure plate
L_launcher = 0.3500 # Length of the launcher

# x-locations:
x_buf =-0.2100 	# was 0.200 #Location of buffer
x_launcher = -14.522 #Location of of launcher start
#x_reservoir_downstream = x_launcher - 2*L_lead_in - L_launcher  # Location of downstream end of reservoir
#x_reservoir_upstream = x_reservoir_downstream - 8.49  # Location of upstream end of reservoir
x_pd =  0.000    # Location of primary diaphragm
x_sd =  9.772   # Location of secondary diaphragm
x_td = 22.237  # (=23.297) Location of tertiary diaphragm
#x_nozzle_inlet = 38.599
#x_tube_end = 40.0  # (=60.379) Location of acceleration tube exit

# nozzle x_coords:
"""
x_n1 = 0 + x_nozzle_inlet
x_n2 = 0.1 + x_nozzle_inlet
x_n3 = 0.2 + x_nozzle_inlet
x_n4 = 0.3 + x_nozzle_inlet
x_n5 = 0.4 + x_nozzle_inlet
x_n6 = 0.5 + x_nozzle_inlet
x_n7 = 0.6 + x_nozzle_inlet
x_n8 = 0.7 + x_nozzle_inlet
x_n9 = 0.8 + x_nozzle_inlet
x_n10 = 0.9 + x_nozzle_inlet
x_n11 = 1.0 + x_nozzle_inlet
x_n12 = 1.1 + x_nozzle_inlet
x_n13 = 1.2 + x_nozzle_inlet
x_n14 = 1.3 + x_nozzle_inlet
x_n15 = 1.4 + x_nozzle_inlet
x_n16 = 1.5 + x_nozzle_inlet
x_n17 = 1.6 + x_nozzle_inlet
x_n18 = 1.7 + x_nozzle_inlet
x_n19 = 1.8 + x_nozzle_inlet
x_n20 = 1.9 + x_nozzle_inlet
x_n21 = 2.0 + x_nozzle_inlet
x_n22 = 2.1 + x_nozzle_inlet
x_n23 = 2.2 + x_nozzle_inlet
x_n24 = 2.3 + x_nozzle_inlet

# nozzle y_coords:
y_n1 = 0.0913*2
y_n2 = 0.0945938*2
y_n3 = 0.1013148*2
y_n4 = 0.1098183*2
y_n5 = 0.1193289*2
y_n6 = 0.1304432*2
y_n7 = 0.1398141*2
y_n8 = 0.1503710*2
y_n9 = 0.1609940*2
y_n10 = 0.1716389*2
y_n11 = 0.1822949*2
y_n12 = 0.1929743*2
y_n13 = 0.2037036*2
y_n14 = 0.2145153*2
y_n15 = 0.2254403*2
y_n16 = 0.2365000*2
y_n17 = 0.2476992*2
y_n18 = 0.2590185*2
y_n19 = 0.2704076*2
y_n20 = 0.2817783*2
y_n21 = 0.2929979*2
y_n22 = 0.3038823*2
y_n23 = 0.3141898*2
y_n24 = 0.3236138*2
"""
# New reservoir arrangement features
# diameters
D_launcher = 0.231
D_launcher_edge = 0.24
D_launcher_bump = 0.265
D_launcher_tube = 0.250
D_launcher_tube_holes = 0.250
D_body_cyl_top = 0.210
D_body_cyl_bot = 0.330
D_block_holes = 0.200
D_block_cavity = 0.4 # changed from 0.460
D_block_new_res_in = 0.225
D_new_res_inlet = 1.2*0.172 # changed from 0.172
D_new_res = 0.68 # changed from 0.71
D_old_res = 0.173
#positions starting from 0 at launcher (-14.552)
l_launcher_holes = 0.050
l_launcher = 0.250
l_launcher_bump = 0.050
l_launcher_tube = 1.100
l_body_cyl_str = 0.164
l_body_cyl_curve = .085
l_block_holes = 0.060
l_block_cavity = 0.171
l_block_new_res = 0.125
l_new_res_inlet = 0.306
l_new_res = 0.720
l_old_res = 11.364

l_list = [l_launcher_holes, l_launcher, l_launcher_bump, l_launcher_tube, l_body_cyl_str, l_body_cyl_curve, l_block_holes, l_block_cavity, l_block_new_res, l_new_res_inlet, l_new_res, l_old_res]

# Diameters:
D_driver        = 0.5000    # Driver tube
D_buffers	= 0.4115    #NEW, EQUIVALENT VOLUME WITH 12 x d=82mm BUFFERS INSTALLED
D_pp            = 0.3000    # Pressure plate
D_orifice	= 0.1410    # Driver orifice plate
D_first_driven  = 0.2000    # 1st driven tube
D_second_driven = 0.1800    # 2nd driven tube
D_accel_tube    = 0.0600    # Acceleration tube
D_bufferplate   = 0.3527    #BUFFER PLATE EQUIVALENT DIAMETER, ACCORDING TO REMAINING AREA

#Define break points of launcher/reservoir - will fix later

add_break_point(x_launcher-np.sum(l_list[0:12]), D_old_res)
add_break_point(x_launcher-np.sum(l_list[0:11])-0.5-0.5, D_old_res)
add_break_point(x_launcher-np.sum(l_list[0:11])+0.2-0.5, D_new_res)
add_break_point(x_launcher-np.sum(l_list[0:10])-0.3-0.5, D_new_res)
add_break_point(x_launcher-np.sum(l_list[0:10])+0.2-0.5, D_old_res)
add_break_point(x_launcher-np.sum(l_list[0:10]), D_old_res)
add_break_point(x_launcher-np.sum(l_list[0:10])+0.2, D_new_res_inlet) # This is the bit we will keep the same from here down
add_break_point(x_launcher-np.sum(l_list[0:9])-0.05, D_new_res_inlet)
#add_break_point(x_launcher-np.sum(l_list[0:9])+0.05, D_block_new_res_in)
#add_break_point(x_launcher-np.sum(l_list[0:8])-0.05, D_block_new_res_in)
add_break_point(x_launcher-np.sum(l_list[0:8])+0.05, D_block_cavity)
add_break_point(x_launcher-np.sum(l_list[0:7])-0.05, D_block_cavity)
add_break_point(x_launcher-np.sum(l_list[0:7])+0.05, D_block_holes)
add_break_point(x_launcher-np.sum(l_list[0:6]), D_block_holes)
#add_break_point(x_launcher-np.sum(l_list[0:6])+0.05, D_body_cyl_top)
#add_break_point(x_launcher-np.sum(l_list[0:5]), D_body_cyl_top)
#add_break_point(x_launcher-np.sum(l_list[0:4]), D_body_cyl_top)
add_break_point(x_launcher-np.sum(l_list[0:4])+0.15, D_launcher_tube)
#add_break_point(x_launcher-np.sum(l_list[0:3]), D_launcher_tube)
#add_break_point(x_launcher-np.sum(l_list[0:3])+0.05, D_launcher_tube)
add_break_point(x_launcher-np.sum(l_list[0:2]), D_launcher_tube)
add_break_point(x_launcher-np.sum(l_list[0:2])+0.05, D_launcher_edge)
add_break_point(x_launcher-l_list[0], D_launcher_edge)
add_break_point(x_launcher-l_list[0]+0.025, D_launcher)
add_break_point(x_launcher, D_launcher)

# Define the tube inside wall profile
add_break_point(x_launcher+0.010, D_driver)
add_break_point(x_buf-0.020, D_driver)
add_break_point(x_buf+0.020, D_buffers)
if False:
    # 2020-06-02 PJ finds the following transtions to and from
    # the buffer-plate diameter are too sharp for L1d4.
    add_break_point(x_buf+0.120-0.020, D_buffers)
    add_break_point(x_buf+0.120+0.020, D_bufferplate)
    add_break_point(x_buf+0.120+0.040-0.020, D_bufferplate)
    add_break_point(x_buf+0.120+0.040+0.020, D_buffers)
add_break_point(x_pd-0.110-0.020, D_buffers)
add_break_point(x_pd-0.110+0.020, D_pp)
add_break_point(x_pd-0.020, D_pp)
add_break_point(x_pd+0.020, D_orifice)
add_break_point(x_pd+0.200-0.020, D_orifice)
add_break_point(x_pd+0.200+0.020, D_first_driven)
add_break_point(x_sd, D_first_driven)
add_break_point(x_sd+L_lead_in, D_second_driven)
add_break_point(x_td, D_second_driven)
add_break_point(x_td+L_lead_in, D_accel_tube)
# add_break_point(x_td+L_lead_in+0.5, D_accel_tube)
# add_break_point(x_nozzle_inlet, D_second_driven)
# add_break_point(x_nozzle_inlet+L_lead_in, D_accel_tube)

"""
add_break_point(x_n1, y_n1)
add_break_point(x_n2, y_n2)
add_break_point(x_n3, y_n3)
add_break_point(x_n4, y_n4)
add_break_point(x_n5, y_n5)
add_break_point(x_n6, y_n6)
add_break_point(x_n7, y_n7)
add_break_point(x_n8, y_n8)
add_break_point(x_n9, y_n9)
add_break_point(x_n10, y_n10)
add_break_point(x_n11, y_n11)
add_break_point(x_n12, y_n12)
add_break_point(x_n13, y_n13)
add_break_point(x_n14, y_n14)
add_break_point(x_n15, y_n15)
add_break_point(x_n16, y_n16)
add_break_point(x_n17, y_n17)
add_break_point(x_n18, y_n18)
add_break_point(x_n19, y_n19)
add_break_point(x_n20, y_n20)
add_break_point(x_n21, y_n21)
add_break_point(x_n22, y_n22)
add_break_point(x_n23, y_n23)
add_break_point(x_n24, y_n24)
"""
