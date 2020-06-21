# heg-495.py
# PJ, 2020-06-13
# Reconstructed from details in report
# Gas-dynamic modelling of the HEG shock tunnel
# PA Jacobs and AD Gardner, 2003.
config.title = "HEG Shock Tunnel, shot 495, Condition V"
#
gm_ideal_air = add_gas_model('ideal-air-gas-model.lua')
gm_he_ar = add_gas_model('he-ar-gas-model.lua') # thermally-perfect but no reactions
massf_he_ar = config.gmodels[gm_he_ar].molef2massf({'He':0.85, 'Ar':0.15})
gm_eq_air = add_gas_model('cea-lut-air.lua')
#
# Tube wall definition follows.
#
# Note that, where ever a restriction is located, we need to have a significant
# length of tube for the gas cells to actually see the correct tube area.
# Also, the transitions between constant-area sections need to be smooth
# enough for the gas-dynamic simulation to be well behaved.
tube.n = 11000 # number of linear segments used for fast interpolation
add_break_point(-3.043, 1.4506)  # upstream-end of new compressed-air reservoir
add_break_point(-0.080, 1.4506)  # start contraction to the holes in the CT
add_break_point(-0.030, 0.6240)  # end contraction, start holes into CT
add_break_point(-0.010, 0.6240)  # end holes into CT
add_break_point( 0.000, 0.5500)  # start of CT proper
add_break_point(33.850, 0.5500)  # start transition to inside of brake cylinder
add_break_point(33.880, 0.4940)  # the inside of the brake cylinder (a.k.a. buffer at UQ)
add_break_point(34.070, 0.4940)  # downstream-end of brake cylinder
add_break_point(34.130, 0.1500)  # beginning of shock tube diameter
add_break_point(34.180, 0.1500)  # start of contraction to ST throat
add_break_point(34.230, 0.1500)  # throat at upstream-end of shock tube (a.k.a. orifice plate)
add_break_point(34.280, 0.1500)  # expanded back to ST diameter
add_break_point(51.240, 0.1500)  # start of contraction to ST throat
add_break_point(51.280, 0.0500)  # entry region to nozzle throat
add_break_point(51.295, 0.0446)  #
add_break_point(51.323, 0.0270)  #
add_break_point(51.343, 0.0220)  # nozzle throat
add_break_point(51.383, 0.0275)  #
add_break_point(51.880, 0.1000)  # arbitrary expansion, not representative of nozzle profile
#
# for loss regions...
add_loss_region(-0.050,  0.000, 0.25)  # entry to compression tube (K_a)
add_loss_region(33.880, 33.980, 0.25)  # start of piston brake (K_b)
add_loss_region(34.130, 34.180, 0.00)  # contraction near main diaphragm (K_c)
add_loss_region(34.180, 34.280, 0.25)  # orifice plate at start of shock tube (K_34)
add_loss_region(38.400, 38.600, 1.50)  # loss region downstream of orifice (K_38)
add_loss_region(51.280, 51.383, 0.25)  # nozzle throat region (K_d)
#
# tube wall temperatures
tube.T_nominal = 292.0
# add_T_patch(26.0, 27.0, 296.0)  # for example
#
# Components that will form the internals of the machine...
#
piston = Piston(mass=275.2, diam=0.550, xL0=0.0, xR0=0.800, vel0=0.0,
                front_seal_f=0.2, front_seal_area=0.020*0.550*math.pi,
                with_brakes=False, x_buffer=33.46)
diaph1 = Diaphragm(x0=34.08, p_burst=50.0e6)
diaph2 = Diaphragm(x0=51.295, p_burst=200.0e3)
#
# Create the gas-path.
left_wall=VelocityEnd(x0=-3.043, vel=0.0)
compressed_air = GasSlug(p=5.36e6, vel=0.0, T=309.0, gmodel_id=gm_ideal_air,
                         ncells=200, to_end_L=False, to_end_R=False, cluster_strength=1.1,
                         viscous_effects=1)
driver_gas = GasSlug(p=72.7e3, vel=0.0, T=292.0, gmodel_id=gm_he_ar, massf=massf_he_ar,
                     ncells=300, to_end_L=False, to_end_R=True, cluster_strength=1.05,
                     viscous_effects=1)
driven_gas = GasSlug(p=70.0e3, vel=0.0, T=292.0, gmodel_id=gm_eq_air,
                     ncells=300, viscous_effects=1)
air_test_section = GasSlug(p=40.0e3, vel=0.0, T=292.0, gmodel_id=gm_eq_air,
                           ncells=50, viscous_effects=1)
right_wall = FreeEnd(x0=51.88)

assemble_gas_path(left_wall, compressed_air, piston, driver_gas, diaph1,
                  driven_gas, diaph2, air_test_section, right_wall)
#
# Other simulation control parameters...
#
config.max_time = 170.0e-3 # the simulation will stop at this time
config.max_step = 500000   # large enough to allow max_time to be reached
config.dt_init = 0.5e-6    # a small enough time step to be stable at the beginning
add_cfl_value(0.0, 0.4)    # should be small enough to cope with diaphragm rupture, etc
add_dt_plot(0.000, 1.00e-3, 0.050e-3) # Most of the compression process is slow.
add_dt_plot(0.154, 0.03e-3, 1.000e-6) # Record more frequently for the shock waves.
#
# Sensor locations...
# Remember that these history locations will be numbered from 0 in the set
# of history data files.  Thus, the nozzle-supply sensor with be numbered 7.
add_history_loc(33.900)   # downstream-end of the compression tube
add_history_loc(34.300)   # SS1, just downstream of ST throat
add_history_loc(34.683)   # SS2
add_history_loc(38.633)   # SS3
add_history_loc(41.228)   # SS4
add_history_loc(44.563)   # SS5
add_history_loc(47.898)   # SS6
add_history_loc(51.200)   # SS7, downstream-end of shock tube (nozzle-supply region)
