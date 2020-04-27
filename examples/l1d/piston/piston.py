# piston.py
# Dimensions match verification case 2 in the 16th AFMC paper (2007)
# "Development of Casbar: a Two-phase Flow Code for the Interior Ballistics
# Problem", R.J Gollan, I.A. Johnston, B.T. O'Flaherty and P.A. Jacobs

config.title = "Ideal piston and gas slug, 2020-04-08"
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

# Set some time-stepping parameters
config.dt_init = 1.0e-6
config.max_time = 40.0e-3
config.max_step = 5000
add_dt_plot(0.0, 1.0e-3, 0.1e-3)
add_history_loc(-8.99) # Upstream-end of driver.
