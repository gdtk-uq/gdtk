# piston-with-valve.py
# Dimensions match verification case 2 in the 16th AFMC paper (2007)
# "Development of Casbar: a Two-phase Flow Code for the Interior Ballistics
# Problem", R.J Gollan, I.A. Johnston, B.T. O'Flaherty and P.A. Jacobs
# Note that the valve and its opening time is arbitrary and is just to
# demonstrate a change in behaviour from the ideal case. PJ 2020-08-15

config.title = "Ideal piston and gas slug plus valve, 2020-08-14"
my_gm = add_gas_model('ideal-air-gas-model.lua')

# Define the tube walls.
add_break_point(-9.0, 0.4)
add_break_point( 7.0, 0.4)
tube.n = 100  # not too many pieces of tube wall

# Create the gas-path.
# Although the gdata name is special, the other variables
# created below are determined by the user for their convenience.
left_end = VelocityEnd(x0=-9.0, vel=0.0)
n_gas_cells = 100
driver_gas = GasSlug(p=100.0e3, vel=0.0, T=348.43, gmodel_id=my_gm,
                     ncells=n_gas_cells, to_end_R=1, cluster_strength=1.1,
                     hcells=n_gas_cells-1)
piston = Piston(mass=1.0, diam=0.4, xL0=-0.005, xR0=0.005, vel0=0.0)
assemble_gas_path(left_end, driver_gas, piston)

v = Valve(x=-1.0, times=[0.020, 0.030], fopen=[0.0, 1.0])

# Set some time-stepping parameters
config.dt_init = 1.0e-6
add_cfl_value(0.0, 0.1) # make small so that the valve interacton is not so bad
config.max_time = 40.0e-3
config.max_step = 5000
add_dt_plot(0.0, 1.0e-3, 0.1e-3)
add_history_loc(-8.99) # Upstream-end of driver.
