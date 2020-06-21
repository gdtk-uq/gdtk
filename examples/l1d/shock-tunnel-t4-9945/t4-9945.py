# t4-9945.py
# PJ, 2012-10-01 L1d3, 2020-05-22 L1d4
# Sarah Razzaqi's 7.6MJ/kg shot with 3mm diaphragm.
config.title = "T4 Tunnel Simulation for run 9945"
#
gm_ideal_air = add_gas_model('ideal-air-gas-model.lua')
gm_he_ar = add_gas_model('he-ar-gas-model.lua') # thermally-perfect but no reactions
massf_he_ar = config.gmodels[gm_he_ar].molef2massf({'He':0.8, 'Ar':0.2})
gm_eq_air = add_gas_model('cea-lut-air.lua')
#
# Tube wall definition follows.
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
add_break_point(-0.20, 0.214)  # start of equivalent duct modelling slots
add_break_point(-0.10, 0.214)  # end of slots.
add_break_point( 0.00, 0.229)  # upstream-end of compression tube
add_break_point(25.85, 0.229)  # downstream-end of compression tube
add_break_point(26.00, 0.067)  # start of restricted entry to shock tube
add_break_point(26.10, 0.067)  # end of restricted entry to shock tube
add_break_point(26.20, 0.076)  # full bore of shock tube
add_break_point(35.90, 0.076)  # downstream-end of shock tube
add_break_point(36.00, 0.025)  # start of nozzle throat
add_break_point(36.05, 0.025)  # end of nozzle throat
add_break_point(37.00, 0.262)  # exit-plane of nozzle, start of test section
#
# for loss regions...
add_loss_region(-0.20, -0.10,  0.5)
add_loss_region(26.00, 26.10,  0.5)
add_loss_region(28.00, 28.50,  1.0)
add_loss_region(36.95, 36.05, 0.25)
#
# tube wall temperatures
tube.T_nominal = 296.0
# add_T_patch(26.0, 27.0, 296.0)  # for example
#
# Components that will form the internals of the machine...
#
piston = Piston(mass=87.0, diam=0.229, xL0=0.0, xR0=0.470, vel0=0.0,
                front_seal_f=0.2, front_seal_area=0.0135*0.229*math.pi,
                with_brakes=True, x_buffer=25.70)
diaph1 = Diaphragm(x0=26.0, p_burst=36.0e6)
diaph2 = Diaphragm(x0=36.0, p_burst=200.0e3)
#
# Create the gas-path.
left_wall=VelocityEnd(x0=-4.87, vel=0.0)
compressed_air = GasSlug(p=2.4e6, vel=0.0, T=300.0, gmodel_id=gm_ideal_air,
                         ncells=200, to_end_L=False, to_end_R=True, cluster_strength=1.2,
                         viscous_effects=1)
driver_gas = GasSlug(p=40.2e3, vel=0.0, T=296.0, gmodel_id=gm_he_ar, massf=massf_he_ar,
                     ncells=200, to_end_L=False, to_end_R=True, cluster_strength=1.05,
                     viscous_effects=1)
driven_gas = GasSlug(p=50.0e3, vel=0.0, T=296.0, gmodel_id=gm_eq_air,
                     ncells=200, viscous_effects=1)
air_test_section = GasSlug(p=40.0e3, vel=0.0, T=296.0, gmodel_id=gm_eq_air,
                           ncells=25, viscous_effects=1)
right_wall = FreeEnd(x0=37.0)

assemble_gas_path(left_wall, compressed_air, piston, driver_gas, diaph1,
                  driven_gas, diaph2, air_test_section, right_wall)
#
# Other simulation control parameters...
#
config.max_time = 268.0e-3 # the simulation will stop at this time
config.max_step = 250000   # large enough to allow max_time to be reached
config.dt_init = 0.5e-6    # a small enough time step to be stable at the beginning
add_cfl_value(0.0, 0.4)    # should be small enough to cope with diaphragm rupture, etc
add_dt_plot(0.0, 2.0e-3, 1.0e-3)    # Most of the compression process is slow.
add_dt_plot(0.251, 0.05e-3, 1.0e-6) # Record more frequently for the shock waves.
#
# Sensor locations...
# Remember that these history locations will be numbered from 0 in the set
# of history data files.  Thus, the nozzle-supply sensor with be numbered 4.
add_history_loc(25.8)   # downstream-end of the compression tube
add_history_loc(30.0)   # shock-tube station 1
add_history_loc(32.0)   # shock-tube station 2
add_history_loc(34.0)   # shock-tube station 3
add_history_loc(35.85)  # downstream-end of shock tube (nozzle-supply region)
add_history_loc(37.0)   # downstream-end of nozzle
