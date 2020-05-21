# shock-he-n2.py
config.title = 'Shock tube, helium driving reacting nitrogen, 2020-05-21'
gm_he = add_gas_model('ideal-helium-gas-model.lua')
gm_n2 = add_gas_model('nitrogen-2sp.lua', 'chem.lua')
config.reacting = True

# Define the tube walls.
add_break_point(0.0, 0.01)
add_break_point(3.0, 0.01)

# Create the gas-path.
left_wall = VelocityEnd(x0=0.0, vel=0.0)
driver_gas = GasSlug(p=1.35e6, vel=0.0, T=3000.0, gmodel_id=gm_he, ncells=300)
interface = GasInterface(x0=0.5)
driven_gas = GasSlug(p=500.0, vel=0.0, T=700.0, massf={'N2':1.0},
                     gmodel_id=gm_n2, ncells=200)
right_wall = VelocityEnd(x0=1.0, vel=0.0)
assemble_gas_path(left_wall, driver_gas, interface, driven_gas, right_wall)

# Set some time-stepping parameters
config.dt_init = 1.0e-8
config.max_time = 80.0e-6
config.max_step = 9000
add_dt_plot(0.0, 2.0e-6, 1.0e-6)
add_history_loc(0.7)
