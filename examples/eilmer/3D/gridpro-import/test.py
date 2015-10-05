# Author: Rowan J. Gollan
# Date: 16-Aug-2012
# Place: Greenslopes, Queensland, Australia

gpro_grid = 'blk.tmp'
gpro_pty = 'blk.tmp.pty'

gdata.title = "Gridpro import test"
gdata.dimensions = 3

u_inf = 900.0
T_inf = 300.0
p_inf = 500.0

select_gas_model(model='ideal gas', species=['air'])

inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=[T_inf], massf=[1.0])
initial = FlowCondition(p=p_inf/3.0, u=0.0, v=0.0, T=[T_inf], massf=[1.0])

gdata.flux_calc = AUSMDV
gdata.dt = 1.0e-6
gdata.cfl = 0.4

# Read in grid from gridpro
grids = read_gridpro_grid(gpro_grid)

# Now loop to set up blocks
blk_list = []
for ib in range(len(grids)):
    blk_list.append(Block3D(grid=grids[ib], fill_condition=initial))

# identify_block_connections()
apply_gridpro_connectivity('blk.tmp.conn', blk_list) 
print type(inflow)
bc_map = {'SUP_IN': inflow,
          'FIXED_T': 450.0 }

apply_gridpro_bcs(gpro_pty, blk_list, bc_map)


