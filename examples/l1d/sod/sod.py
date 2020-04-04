# sod.py
config.title = 'Sods ideal shock tube, 2020-04-04'
my_gm = add_gas_model('ideal-air-gas-model.lua')

# Define the tube walls.
add_break_point(0.0, 0.01)
add_break_point(3.0, 0.01)

# Create the gas-path.
left_wall = VelocityEnd(x0=0.0, vel=0.0)
driver_gas = GasSlug(p=100.0e3, vel=0.0, T=348.4, gmodel_id=my_gm, nn=200)
interface = GasInterface(x0=0.5)
driven_gas = GasSlug(p=10.0e3, vel=0.0, T=278.7, gmodel_id=my_gm, nn=100)
right_wall = VelocityEnd(x0=1.0, vel=0.0)
assemble_gas_path(left_wall, driver_gas, interface, driven_gas, right_wall)

# Set some time-stepping parameters
config.dt_init = 1.0e-7
config.max_time = 0.6e-3
config.max_step = 5000
add_dt_plot(0.0, 10.0e-6, 5.0e-6)
add_history_loc(0.7)
