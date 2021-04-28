# File: x2_condition1.py
#
# This is an L1d4 simulation of the X2 expansion tube without 
# a secondary driver section for the 2020 Nonequilibrium ARC Discovery
#
# 14 April 2021: Adapted to the first air X2 condition, Carolyn Jacobs
#
# Some history:
# 15 Feb 2021: Simplified by Carolyn Jacobs
# Feb 2021: Adapted by Arlo Leckie (arlo.leckie@uq.net.au) - refinment of 
#           undergraduate thesis simulation
# This input file is based on work done by David Gildfind and Umar Sheikh.
# Based on original simulation by Chris James (c.james4@uq.edu.au) 

##############################################################################
# Preamble information                                                       #
##############################################################################

import math

config.title = 'X2 expansion mode, 2020 nonequilibrium DP condition 1'

# Scale factor to simplify grid resolution study
mesh_scale = 10

# Flag to control the viscous settings (0 = inviscid, 1 = viscous)
vis_flag = 0

# Fill pressures 
p_res_fill      = 6.85e6      # Reservoir fill pressure [Pa]
p_drv_fill      = 92.8e3      # Driver fill pressure [Pa]
prim_dia_burst  = 27.9e6      # Primary diaphragm burst pressure [Pa]
sec_dia_burst   = 18.0e3      # Secondary diaphragm burst pressure [Pa]
shk_tube_fill   = 3000.0      # Shock tube fill pressure[Pa]
acc_tube_fill   = 10.0	      # Acceleration tube fill pressure [Pa]

# Secondary diaphragm burst pressure is from the burst pressure of a single 
# aluminium foil diaphragm.

# Fixed geometry information
x_res_start = -3.890                    # x-loc of the reservoir start [m]
x_buffer = 4.600                        # x-loc of buffer estimate [m]
x_pri_diaphragm = 4.810                 # x-loc of the primary diaphragm [m]
x_sec_diaphragm = x_pri_diaphragm+3.418 # x-loc of secondary diaphragm [m]

# Fixed transducer locations
x_sd1 = x_pri_diaphragm+2.577
x_sd2 = x_pri_diaphragm+2.810
x_sd3 = x_pri_diaphragm+3.043
x_st1 = x_pri_diaphragm+4.231
x_st2 = x_pri_diaphragm+4.746
x_st3 = x_pri_diaphragm+5.260
x_at1 = x_pri_diaphragm+6.437
x_at2 = x_pri_diaphragm+6.615
x_at3 = x_pri_diaphragm+6.796
x_at4 = x_pri_diaphragm+7.590
x_at5 = x_pri_diaphragm+7.846
x_at6 = x_pri_diaphragm+8.096

# Nominal end of straight section (adjusted later depending on configuration)
x_end_acceleration = x_pri_diaphragm+8.585

##############################################################################
# Tube geometry                                                              #
##############################################################################

# Reservoir and compression tube:
add_break_point(x_res_start, 0.3160)
add_break_point(-0.990, 0.3160)
add_break_point(-0.970, 0.2440)
add_break_point(-0.370, 0.2440)
add_break_point(-0.350, 0.1600)
add_break_point(-0.157, 0.1600)
add_break_point(-0.010, 0.2568)

# Tunnel downstream of compression tube:
add_break_point(4.600, 0.2568) # Beginning of area change at PD
add_break_point(4.700, 0.0850) # End of area change to PD
# Orifice plate is placed here for the other driver configurations, just after 
# the PD which is at 4.810

# EITHER Add the nozzle: (refer Michael Scott PhD Thesis, 2006, Table 5.3)
# THIS IS THE AT EXIT, SET ABOVE. ADD OFFSET OF 0.958 FROM HERE
noz_offset = -0.010
nozzle_start = x_end_acceleration+noz_offset
nozzle_exit = nozzle_start+1.4
add_break_point(nozzle_start, 0.085) 
add_break_point(nozzle_start+0.023422, 0.08696)
add_break_point(nozzle_start+0.052994, 0.089496)
add_break_point(nozzle_start+0.089767, 0.092956)
add_break_point(nozzle_start+0.134836, 0.097684)
add_break_point(nozzle_start+0.189295, 0.103878)
add_break_point(nozzle_start+0.254202, 0.111644)
add_break_point(nozzle_start+0.330529, 0.121042)
add_break_point(nozzle_start+0.419119, 0.132066)
add_break_point(nozzle_start+0.520628, 0.144532)
add_break_point(nozzle_start+0.635471, 0.157942)
add_break_point(nozzle_start+0.763763, 0.17141)
add_break_point(nozzle_start+0.905261, 0.183712)
add_break_point(nozzle_start+1.059306, 0.193498)
add_break_point(nozzle_start+1.247690, 0.199654)
add_break_point(nozzle_exit, 0.20168)
#
# OR Add a straight adaptor to the dumptank:
# add_break_point(13.789, 0.0850, 1) # Analysis tube exit. 
#
##############################################################################
# Notes about the tube configuration:                                        #
#                                                                            #
# Configuration may have mylar diaphragms at x = 8.228m, x = 10.786	         #
# i.e. these are the capstan locations where the normally used tubes are     #
#      joined.                                                               #
#                                                                            #
# The normal fixed length 85mm diameter tubes end at x=13.169m.              #
#                                                                            #
# If the facility is to be run in tube mode, then a 0.620m long dumptank     #
# adaptor is installed, in which case the tube ends at x=13.169+0.620=13.789m#
#                                                                            #
# If the facility is to be run with the nozzle, 0.226m of the nozzle is      #
# straight, therefore the straight section ends at 13.169+0.226=13.395m.     #    
# The contoured part of the nozzle is 1.4m long, therefore the nozzle exit   #
# is located at x=13.395+1.4=14.795m.                                        #
##############################################################################

##############################################################################
# Gas model details                                                          #
##############################################################################

gm_ideal_air = add_gas_model('ideal-air-gas-model.lua')
gm_he_ar = add_gas_model('he-ar-gas-model.lua')
massf_he_ar = config.gmodels[gm_he_ar].molef2massf({'He':0.8, 'Ar':0.2})
gm_eq_air = add_gas_model('cea-lut-air.lua')

##############################################################################
# Construct the gas path                                                     #
##############################################################################

# Ambient air temperature
T_amb = 22.0+273.15

# Create the gas-path
left_wall = VelocityEnd(x0=x_res_start, 
                        vel=0.0)
res_gas = GasSlug(gmodel_id=gm_ideal_air, 
                  p=p_res_fill, 
                  T=T_amb, 
                  vel=0.0,
                  ncells=10*mesh_scale, 
                  viscous_effects=vis_flag, 
                  adiabatic=False,
                  hcells=1,
                  label='compressed air in reservoir')
piston = Piston(mass=10.524, 
                diam=0.2568, 
                xL0=0.0, 
                xR0=0.221, 
                vel0=0.0,
                front_seal_f=0.2, 
                front_seal_area=0.020*0.2568*math.pi,
                x_buffer=4.5895, 
                on_buffer=0,
                label='lightweight piston')
driver_gas = GasSlug(gmodel_id=gm_he_ar, 
                     p=p_drv_fill, 
                     T=T_amb, 
                     vel=0.0,
                     massf=massf_he_ar,
                     ncells=10*mesh_scale, 
                     cluster_strength=1.05,
                     to_end_R=True,
                     viscous_effects=vis_flag, 
                     adiabatic=False,
                     hcells=1,
                     label='80/20 by mole He-Ar driver gas')
primary_diaphragm = Diaphragm(x0=x_pri_diaphragm, 
                              p_burst=prim_dia_burst, 
                              state=0)
test_gas = GasSlug(gmodel_id=gm_eq_air, 
                   p=shk_tube_fill, 
                   T=T_amb, 
                   vel=0.0, 
                   ncells=10*mesh_scale, 
                   cluster_strength=1.05,
                   to_end_R=True,
                   viscous_effects=vis_flag, 
                   adiabatic=False,
                   hcells=1,
                   label='test gas')
secondary_diaphragm = Diaphragm(x0=x_sec_diaphragm, 
                                p_burst=sec_dia_burst, 
                                state=0)
accelerator_gas = GasSlug(gmodel_id=gm_eq_air, 
                          p=acc_tube_fill, 
                          T=T_amb, 
                          vel=0.0,
                          ncells=10*mesh_scale, 
                          cluster_strength=1.05,
                          to_end_L=True,
                          viscous_effects=vis_flag, 
                          adiabatic=False,
                          hcells=1,
                          label='accelerator gas')
right_free = FreeEnd(x0=nozzle_exit)

# Assemble the gas path
assemble_gas_path(left_wall, res_gas, piston, driver_gas, primary_diaphragm, 
                  test_gas, secondary_diaphragm, accelerator_gas, right_free)

# Add loss regions
L_lead_in = 0.1
add_loss_region(0.0,L_lead_in,0.1)              # Piston launch tube losses
add_loss_region(x_buffer,x_pri_diaphragm,0.7)   # Primary diaphragm losses
add_loss_region(x_sec_diaphragm,x_sec_diaphragm+L_lead_in,0.2)

##############################################################################
# General solution parameters                                                #
##############################################################################

# Time-stepping parameters we might need to tweak
t_finish=0.030    # Simulation end time [s]
t_switch=0.023    # Simulation switch to finer time steps [s]
t_fine=1.0e-5*5.0 # Size of finer time steps [s]

# General time-stepping parameters
config.dt_init = 1.0e-10
config.max_time = t_finish
config.max_step = 25000000
add_cfl_value(0.0, 0.25)
add_dt_plot(0.0, 2.0e-4, 2.0e-4)
add_dt_plot(t_switch, t_fine, t_fine/10)

# Define history locations
add_history_loc(4.600)              # Compression tube pressure immediately upstream of primary diaphragm. 0
add_history_loc(x_pri_diaphragm)    # Compression tube pressure at primary diaphragm. 1
add_history_loc(x_sd1)              # PCB transducer sd1. 2
add_history_loc(x_sd2)              # PCB transducer sd2. 3
add_history_loc(x_sd3)              # PCB transducer sd3. 4
add_history_loc(x_st1)              # PCB transducer st1. 5
add_history_loc(x_st2)              # PCB transducer st2. 6
add_history_loc(x_st3)              # PCB transducer st3. 7
add_history_loc(x_at1)              # PCB transducer at1. 8
add_history_loc(x_at2)              # PCB transducer at2. 9
add_history_loc(x_at3)              # PCB transducer at3. 10
add_history_loc(x_at4)              # PCB transducer at4. 11
add_history_loc(x_at5)              # PCB transducer at5. 12
add_history_loc(x_at6)              # PCB transducer at6. 13
add_history_loc(nozzle_start)       # End of straight tube part of facility. 14
add_history_loc(nozzle_exit)        # Nozzle exit 15
