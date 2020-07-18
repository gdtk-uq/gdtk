# t4-11315.py
# Ryan Whitside, Jen Kunz, Wilson Chan and Peter J.
# This l1d4 input script built by PJ, 2020-07-17, to match
# the 2019-03-21 L1d3 simulations by Wilson and Ryan.
config.title = "T4 Tunnel Simulation for run 11315, nitrogen driver"
#
gm_ideal_air = add_gas_model('ideal-air-gas-model.lua')
gm_ideal_nitrogen = add_gas_model('ideal-nitrogen-gas-model.lua')
gm_eq_air = add_gas_model('cea-lut-air.lua')
#
# Tube wall geometry from accompanying file.
exec(open('./t4-tube.py').read())
# Loss resions are tuned empirically.
add_loss_region(-0.20, -0.10, 25.0)  # launcher slots
add_loss_region(26.00, 26.10, 0.8)   # primary diaphragm station
                                     # (assumed to be re-entrant inlet type)
add_loss_region(26.10, 36.06, 0.35)  # shock tube losses
add_loss_region(36.06, 36.10, 0.50)  # nozzle contraction and throat
#
# tube wall temperatures
tube.T_nominal = 296.0
# add_T_patch(26.0, 27.0, 296.0)  # for example
#
# Components that will form the internals of the machine...
#
# The standard T4 piston (weight known at 2019-03-21).
# Note that even though the brakes were physically installed,
# they were oiled and hence assumed to be non-functional.
piston = Piston(mass=90.05, diam=0.229, xL0=0.0, xR0=0.470, vel0=0.0,
                front_seal_f=0.2, front_seal_area=0.0137*0.229*math.pi,
                with_brakes=True, brakes_friction_force=0.0, x_buffer=25.87)
# Primary (steel) diaphragm.
diaph1 = Diaphragm(x0=26.0, p_burst=25.2e6, dxL=0.1, dxR=0.1)
# Secondary (Mylar) diaphragm - expected to burst instantaneously
# when shock reflects off the end of the shock tube.
diaph2 = Diaphragm(x0=36.08, p_burst=1.0e6)
#
# Create the gas-path.
left_wall = VelocityEnd(x0=-4.87, vel=0.0)
compressed_air = GasSlug(p=2.40e6, vel=0.0, T=296.0, gmodel_id=gm_ideal_air,
                         ncells=500, to_end_L=False, to_end_R=True, cluster_strength=1.2,
                         viscous_effects=1)
driver_gas = GasSlug(p=81.8e3, vel=0.0, T=296.0, gmodel_id=gm_ideal_nitrogen,
                     ncells=500, to_end_L=False, to_end_R=True, cluster_strength=1.05,
                     viscous_effects=1)
driven_gas = GasSlug(p=270.0e3, vel=0.0, T=296.0, gmodel_id=gm_eq_air,
                     ncells=500, viscous_effects=1)
air_test_section = GasSlug(p=40.0e3, vel=0.0, T=296.0, gmodel_id=gm_eq_air,
                           ncells=50, viscous_effects=1)
right_wall = FreeEnd(x0=37.0)

assemble_gas_path(left_wall, compressed_air, piston, driver_gas, diaph1,
                  driven_gas, diaph2, air_test_section, right_wall)
#
# Other simulation control parameters...
#
config.max_time = 328.0e-3 # the simulation will stop at this time
config.max_step = 50000000 # large enough to allow max_time to be reached
config.dt_init = 0.5e-6    # a small enough time step to be stable at the beginning
add_cfl_value(0.0, 0.4)    # should be small enough to cope with diaphragm rupture, etc
add_dt_plot(0.000, 2.00e-3, 1.0e-3) # Most of the compression process is slow.
add_dt_plot(0.250, 0.05e-3, 1.0e-6) # Record more frequently for the shock waves.
