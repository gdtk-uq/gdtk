-- cht-cyl.lua
--
-- Conjugate heat transfer valdiation test case taken from:
--    Time-adaptive, loosely coupled strategy for conjugate heat transfer problems in hypersonic flows
--    Zhang et al. (2014)
--    Journal of Thermophysics and Heat Transfer
--
-- Kyle D. (April 2020)
-- 

config.title = "Mach 6.47 air flow over a hollow cylinder."
config.axisymmetric = false

-- free stream conditions
nsp, nmodes, gmodel = setGasModel('ideal-air-gas-model.lua')
U_inf = 2015.43 -- m/s
T_inf = 241.5 -- K
Twall = 294.4 -- K
ReL = 1.31e06 -- m^-1
R = 287.0 -- kJ/kg.K
gas = GasState:new{gmodel}
gas.T = T_inf
gmodel:updateTransCoeffs(gas)
p_inf = (ReL*gas.mu)/U_inf * (R*T_inf)
inflow = FlowState:new{p=p_inf, T=T_inf, velx=U_inf}
initial = inflow

-- Set up the geometry for defining the grid
Ro = 3.81e-02 -- outer nose radis, metres
Ri = 2.52e-02 -- inner nose radius, metres

-- define FluidBlocks
n0i = 100
n0j = 50

a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=-Ro, y=0.0}
c = Vector3:new{x=0.0, y=Ro}
d = { Vector3:new{x=-1.5*Ro, y=0.0}, Vector3:new{x=-1.5*Ro, y=Ro},
      Vector3:new{x=-Ro, y=2*Ro}, Vector3:new{x=0.0, y=3*Ro} }

sphere_edge = Arc:new{p0=b, p1=c, centre=a}
psurf0 = makePatch{north=Line:new{p0=d[#d], p1=c}, south=Line:new{p0=d[1], p1=b},
		  east=sphere_edge, west=Bezier:new{points=d}}
cf_radial0 = RobertsFunction:new{end0=false, end1=true, beta=1.2}

grid0 = StructuredGrid:new{psurface=psurf0, niv=n0i+1, njv=n0j+1,
			  cfList={north=cf_radial0, south=cf_radial0}}

blk0 = FBArray:new{grid=grid0, initialState=initial, label='fluidblk',
                       bcList={north=OutFlowBC_Simple:new{},
                               east=None,
                               south=WallBC_WithSlip:new{},
                               west=InFlowBC_Supersonic:new{flowState=inflow}},
                       nib=2, njb=2}

-- define SolidBlocks
n1i = 25

e = Vector3:new{x=-Ri, y=0.0}
f = Vector3:new{x=0.0, y=Ri}

sphere_outer_edge = Arc:new{p0=b, p1=c, centre=a}
sphere_inner_edge = Arc:new{p0=e, p1=f, centre=a}
psurf1 = makePatch{north=Line:new{p0=c, p1=f}, south=Line:new{p0=b, p1=e},
		  east=sphere_inner_edge, west=sphere_outer_edge}
cf_radial1 = RobertsFunction:new{end0=true, end1=false, beta=1.1}

grid1 = StructuredGrid:new{psurface=psurf1, niv=n1i+1, njv=n0j+1,
                           cfList={north=cf_radial1, south=cf_radial1}}

blk1 = SolidBlockArray{grid=grid1, initTemperature=Twall, label='solidblk',
                       properties={rho=8030.0, k=16.24, Cp=502.48},
                       bcList={north=SolidAdiabaticBC:new{},
                               east=SolidFixedTBC:new{T=Twall},
                               south=SolidConstantFluxBC:new{},
                               west=None},
                       nib=2, njb=2}

identifyBlockConnections()

-- every MPI task needs a FluidBlock, so make sure you distribute the blocks accordingly
-- e.g. here we have 4 FluidBlocks + 4 SolidBlocks divided amongst 4 MPI tasks
mpiTasks = mpiDistributeBlocks{ntasks=4, dist="load-balance", preassign={[0]=1}}

-- Now set some configuration options
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.apply_limiter = true
config.extrema_clipping = true

config.viscous = true
config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"
config.viscous_signal_factor = 0.1

config.gasdynamic_update_scheme = 'pc'
config.cfl_value = 0.4
config.dt_init = 1.0e-12
config.max_time = 0.001 -- seconds
config.max_step = 8e7
config.dt_plot = config.max_time/100.0

-- loads settings (adjacent_to_solid is a special internal boundary name for coupled fluid/solid boundaries)
config.boundary_groups_for_loads = "adjacent_to_solid"
config.write_loads = true
config.dt_loads = config.max_time/500
