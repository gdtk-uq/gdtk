# he-air.py
config.title = 'OU shock tunnel, helium driving ideal air, experiment D'
gm_he = add_gas_model('ideal-helium.gas')
gm_air = add_gas_model('ideal-air.gas')

# Tunnel-wall diameters.
add_break_point(-2.400, 0.400) # upstream-end of driver
add_break_point(-0.500, 0.400) # downstream-end of driver
add_break_point(-0.400, 0.065) # entry to valve section
add_break_point(-0.100, 0.065) # valve throat
add_break_point(0.000, 0.100) # upstream-end of shock tube
add_break_point(7.000, 0.100) # downstream-end of shock tube
add_break_point(7.050, 0.0345) # nozzle throat start
add_break_point(7.100, 0.0345) # nozzle throat end
# The expansion through the nozzle is poorly modelled; use Eilmer for that.
add_break_point(7.950, 0.197) # exit of nozzle

add_loss_region(-0.400, 0.100, 1.0) # at valve body
add_loss_region( 7.000, 7.100, 0.5) # at nozzle throat
tube.T_nominal = 293.0
add_vf_patch(7.100, 7.950, 0.0) # no viscous effects in nozzle expansion

# Create the gas-path.
left_end = VelocityEnd(x0=-2.400, vel=0.0)
driver_gas = GasSlug(p=2.100e6, vel=0.0, T=293.0, gmodel_id=gm_he, ncells=1200,
                     viscous_effects=1, to_end_R=True, cluster_strength=1.1)
contact = GasInterface(x0=0.0)
driven_gas = GasSlug(p=19.1e3, vel=0.0, T=293.0, gmodel_id=gm_air, ncells=1200,
                     viscous_effects=1)
diaph = Diaphragm(x0=7.000, p_burst=150.0e3)
# Note that the dump_tank_gas is initially at fairly-high pressure.
# This is to mitigate the cells being strongly shocked and the simulation
# running very slowly.  We otherwise ignore the state of this gas.
dump_tank_gas = GasSlug(p=4.0e3, T=293.0, ncells=60, gmodel_id=gm_air,
                        viscous_effects=1, to_end_L=True, cluster_strength=1.1)
right_end = FreeEnd(x0=7.950)
assemble_gas_path(left_end, driver_gas, contact, driven_gas,
                  diaph, dump_tank_gas, right_end)
v = Valve(x=-0.200, times=[0.000, 0.002], fopen=[0.0, 1.0])

# Set some time-stepping parameters
config.dt_init = 0.5e-6
add_cfl_value(0.0, 0.40)
config.max_time = 25.0e-3
config.max_step = 400000
add_dt_plot(0.0, 100.0e-6, 10.0e-6)
add_history_loc(5.124) # p-sensor 667
add_history_loc(6.949) # p-sensor 668
add_history_loc(8.050) # sphere p-sensor 666
