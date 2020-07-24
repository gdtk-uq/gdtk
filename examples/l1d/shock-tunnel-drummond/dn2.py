# dn2.py
config.title = "Drummond tunnel M4 nozzle p4=3.25MPa N2, p1=30kPa N2"
gm_nitrogen = add_gas_model('thermally-perfect-nitrogen.lua')

# Tube wall geometry from accompanying file.
exec(open('./drummond-tube.py').read())
add_loss_region(-3.050, -3.000, 0.5) # at steel-diaphragm station
add_loss_region( 0.050,  0.120, 0.5) # at nozzle throat
tube.T_nominal = 296.0
add_vf_patch(0.1, 0.35, 0.0) # no viscous effects in nozzle expansion

# Create the gas-path.
left_end = VelocityEnd(x0=-3.785, vel=0.0)
driver_gas = GasSlug(p=3.25e6, T=296.0, gmodel_id=gm_nitrogen,
                     ncells=150, viscous_effects=1)
# Note that there is no steel diaphragm in this simulation.
cs = GasInterface(x0=-3.015)
test_gas = GasSlug(p=30.0e3, T=296.0, ncells=300, gmodel_id=gm_nitrogen,
                   viscous_effects=1)
diaph = Diaphragm(x0=0.10, p_burst=150.0e3)
dump_tank_gas = GasSlug(p=4.0e3, T=296.0, ncells=6, gmodel_id=gm_nitrogen,
                        viscous_effects=1, to_end_L=True, cluster_strength=1.1)
right_end = FreeEnd(x0=0.3)
assemble_gas_path(left_end, driver_gas, cs, test_gas,
                  diaph, dump_tank_gas, right_end)

# Set some time-stepping parameters
config.dt_init = 0.5e-6
add_cfl_value(0.0, 0.40)
config.max_time = 8.0e-3
config.max_step = 100000
add_dt_plot(0.0, 30.0e-6, 2.0e-6)
