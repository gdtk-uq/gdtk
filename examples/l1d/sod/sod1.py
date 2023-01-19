# sod1.py
# Demonstrate the use of a variable initial state for a gas slug.

config.title = 'Sods ideal shock tube with one slug, 2023-01-19'
my_gm = add_gas_model('ideal-air-gas-model.lua')

# Define the tube walls.
add_break_point(0.0, 0.01)
add_break_point(3.0, 0.01)

def my_fun(x):
    fs = {'vel':0.0, 'massf':[1.0,]}
    if x < 0.5:
        # High-pressure driver gas
        fs['p'] = 100.0e3
        fs['T'] = 348.4
    else:
        # Low-pressure driven gas
        fs['p'] = 10.0e3
        fs['T'] = 278.7
    return fs

# Create the gas-path.
left_wall = VelocityEnd(x0=0.0, vel=0.0)
all_gas = GasSlug(initial_fs_fun=my_fun, gmodel_id=my_gm, ncells=200)
right_wall = VelocityEnd(x0=1.0, vel=0.0)
assemble_gas_path(left_wall, all_gas, right_wall)

# Set some time-stepping parameters
config.dt_init = 1.0e-7
config.max_time = 0.6e-3
config.max_step = 5000
add_dt_plot(0.0, 10.0e-6, 5.0e-6)
add_history_loc(0.7)
