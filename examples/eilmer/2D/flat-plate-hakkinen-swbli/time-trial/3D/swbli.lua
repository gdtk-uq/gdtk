-- swbli.lua
-- Anand V, 10-October-2015 and Peter J, 2016-11-02
-- Model of Hakkinen et al's 1959 experiment.

config.title = "Shock Wave Boundary Layer Interaction"
print(config.title)
config.dimensions = 3

-- Flow conditions to match those of Figure 6: pf/p0=1.4, Re_shock=2.96e5
p_inf = 6205.0 -- Pa
u_inf = 514.0 -- m/s
T_inf = 164.4 -- degree K

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
inflow = FlowState:new{p=p_inf, velx=u_inf, T=T_inf}

--------------------- Importing GridPro Grid ---------------------
gpro_grid = 'grid.tmp'
gpro_pty = 'grid.tmp.pty'
grids = importGridproGrid(gpro_grid, 1e-3)

--------------------- Block Setup ---------------------
blk_list = {}

for i,g in ipairs(grids) do
   blk_list[#blk_list+1] = FluidBlock:new{grid=g, initialState=inflow}
end

applyGridproConnectivity('grid.tmp.conn', blk_list)

--------------------- Boundary Conditions ---------------------
bc_map = {SUP_IN=InFlowBC_Supersonic:new{flowState=inflow},
	      SLIP_WALL=WallBC_WithSlip:new{},
	      ADIABATIC=WallBC_NoSlip_Adiabatic:new{},
	      EXTRAPOLATE_OUT=OutFlowBC_Simple:new{}}

applyGridproBoundaryConditions(gpro_pty, blk_list, bc_map)

--------------------- MPI Setup ---------------------
mpiDistributeBlocks{ntasks=4, dist="load-balance"}

mm = 1.0e-3 -- metres per mm
L2 = 90.0*mm

config.gasdynamic_update_scheme = "euler"
config.flux_calculator = 'adaptive'
config.viscous = true
config.spatial_deriv_calc = 'divergence'
config.cfl_value = 0.5
config.max_time = 5.0*L2/u_inf -- time in flow lengths
config.max_step = 200000
config.dt_init = 1.0e-8
config.dt_plot = config.max_time/1
