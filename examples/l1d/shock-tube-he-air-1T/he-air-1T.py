# he-air-1T.py
config.title = 'Shock tube, helium driving 1T air, 2023-01-27'
gm_he = add_gas_model('ideal-helium-gas-model.lua')
gm_air = add_gas_model('air-5sp-1T-gas-model.lua', 'air-5sp-1T-chem.lua')
config.reacting = True

# Define the tube walls.
add_break_point(-3.0, 0.01)
add_break_point(9.0, 0.01)

# Create the gas-path.
left_wall = VelocityEnd(x0=-2.0, vel=0.0)
driver_gas = GasSlug(p=0.800e6, vel=0.0, T=800.0, gmodel_id=gm_he, ncells=1200)
interface = GasInterface(x0=0.0)
driven_gas = GasSlug(p=290.0, vel=0.0, T=293.0, massf={'N2':0.78,'O2':0.22},
                     gmodel_id=gm_air, ncells=1200)
right_wall = VelocityEnd(x0=9.0, vel=0.0)
assemble_gas_path(left_wall, driver_gas, interface, driven_gas, right_wall)

# Set some time-stepping parameters
config.dt_init = 1.0e-8
config.max_time = 3.0e-3
config.max_step = 90000
add_dt_plot(0.0, 20.0e-6, 1.0e-6)
add_history_loc(1.0)
add_history_loc(1.1)
add_history_loc(2.0)
add_history_loc(2.1)
add_history_loc(3.0)
add_history_loc(3.1)
add_history_loc(4.0)
add_history_loc(4.1)
add_history_loc(5.0)
add_history_loc(5.1)
add_history_loc(6.0)
add_history_loc(6.1)
add_history_loc(7.0)
add_history_loc(7.1)
add_history_loc(8.0)
add_history_loc(8.1)
