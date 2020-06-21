# x3r-561.py
#
# Simulation of X3R reflected shock tunnel with full piston dynamics.
# Sam Stennett, 2019.
# Adapted to L1d4 by PJ 2020-06-02

import math
import numpy as np

config.title = "X3R Candidate Condition, no rupture of diaphragm."
#
gm_ideal_air = add_gas_model('ideal-air-gas-model.lua')
gm_he_ar = add_gas_model('he-ar-gas-model.lua') # thermally-perfect but no reactions
massf_he_ar = config.gmodels[gm_he_ar].molef2massf({'He':0.5, 'Ar':0.5})
gm_eq_air = add_gas_model('cea-lut-air.lua')
#
# Tube wall geometry from accompanying file.
exec(open('./x3r-tube.py').read())
#
# for loss regions...
add_loss_region(-14.9, -14.872, 0.1) # loss factor over launch tube
# add_loss_region(x_buf, x_pd+L_lead_in, 0.5)  # primary diaphragm
# add_loss_region(x_sd, x_sd+L_lead_in, 0.20)  # secondary diaphragm
#
# tube wall temperatures
T_amb = 25.0 + 273.15
tube.T_nominal = T_amb
#
# Components that will form the internals of the machine...
#
left_wall = VelocityEnd(x0=x_launcher-np.sum(l_list), vel=0.0)
reservoir_gas = GasSlug(p=4.2e6, T=T_amb, gmodel_id=gm_ideal_air,
                        ncells=600, viscous_effects=1,
                        label='compressed air in reservoir')
piston = Piston(mass=280.0, diam=D_driver, xL0=x_launcher, xR0=x_launcher+L_piston,
                x_buffer=x_buf-L_piston/2, vel0=0.0,
                front_seal_f=0.0, front_seal_area=0.0206*D_driver*math.pi,
                back_seal_f=0.1, back_seal_area=0.0206*D_driver*math.pi,
                label="Piston - 280.0 kg version")
driver_one_gas = GasSlug(p=90.0e3, T=T_amb, gmodel_id=gm_he_ar, massf=massf_he_ar,
                         ncells=300, viscous_effects=1,
                         to_end_L=True, to_end_R=True, cluster_strength=1.05,
                         label='50/50 helium-argon driver gas')
primary_diaphragm = Diaphragm(x0=x_pd, p_burst=75.0e6) # Large pressure --> no rupture.
test_gas = GasSlug(p=33.1e3, T=T_amb, gmodel_id=gm_eq_air,
                   ncells=50, # small number because diaphragm will not rupture
                   viscous_effects=0,
                   to_end_L=True, to_end_R=True, cluster_strength=1.02,
                   label='test gas')
right_free = FreeEnd(x0=x_td+L_lead_in)

assemble_gas_path(left_wall, reservoir_gas, piston, driver_one_gas,
                  primary_diaphragm, test_gas, right_free)
#
# Other simulation control parameters...
#
config.max_time = 0.600    # long enough to get two bounces of the piston
config.max_step = 2500000  # the simulation will stop here, dt becomes too small to make progress
config.dt_init = 0.5e-7    # a small enough time step to be stable at the beginning
add_cfl_value(0.0, 0.25)   # should be small enough to cope with diaphragm rupture, etc
# We want a uniform record of the compression process,
# since we are not letting the diaphragm rupture.
add_dt_plot(0.0, 1.0e-3, 0.1e-3) # The compression process is slow.
#
# Specify history locations
add_history_loc(-0.250)       # 0 - ct (compression tube), upstream of buffer
add_history_loc(-0.175)       # 1 - ct (compression tube), within buffer space
# Shock-tube sensors
add_history_loc(1.649)        # 2 - st1
add_history_loc(6.504)        # 3 - st2
add_history_loc(7.955)        # 4 - st3
add_history_loc(9.519)        # 5 - st4
# Acceleration-tube sensors
add_history_loc(16.123)       # 6 - st5
add_history_loc(19.217)       # 7 - st6
add_history_loc(20.715)       # 8 - st7
add_history_loc(22.207)       # 9 - st8
