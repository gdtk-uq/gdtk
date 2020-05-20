# exptube.py
config.title = 'Expansion-tube demo, 2020-05-20'
my_gm = add_gas_model('ideal-air-gas-model.lua')

# Define the tube wall.
add_break_point(0.0, 0.02)
add_break_point(3.0, 0.02)

# Create the gas-path.
left_wall = VelocityEnd(x0=0.0, vel=0.0)
driver_gas = GasSlug(p=1.0e6, vel=0.0, T=3000.0, gmodel_id=my_gm, ncells=200)
interface = GasInterface(x0=0.5)
driven_gas = GasSlug(p=5.0e3, vel=0.0, T=300.0, gmodel_id=my_gm, ncells=100)
dia = Diaphragm(x0=1.0, p_burst=10.0e3)
accel_gas = GasSlug(p=100.0, vel=0.0, T=300.0, gmodel_id=my_gm, ncells=300)
free_end = FreeEnd(x0=2.5)
assemble_gas_path(left_wall, driver_gas, interface, driven_gas,
                  dia, accel_gas, free_end)

# Set some time-stepping parameters
config.dt_init = 1.0e-8
config.max_time = 1.0e-3
config.max_step = 9000
add_dt_plot(0.0, 10.0e-6, 1.0e-6)
add_history_loc(0.90)
add_history_loc(2.45)
