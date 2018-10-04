-- swbli.lua
-- Anand V, 10-October-2015 and Peter J, 2016-11-02
-- Model of Hakkinen et al's 1959 experiment.

config.title = "Shock Wave Boundary Layer Interaction"
print(config.title)
config.dimensions = 2

-- Flow conditions to match those of Figure 6: pf/p0=1.4, Re_shock=2.96e5
p_inf = 6205.0 -- Pa
u_inf = 514.0 -- m/s
T_inf = 164.4 -- degree K

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
inflow = FlowState:new{p=p_inf, velx=u_inf, T=T_inf}

-- Flow domain.
--
--   y
--   ^   a1---b1---c1---d1   Shock generator
--   |   |    |     |    |
--   |   | 0  |  1  | 2  |   patches
--   |   |    |     |    |
--   0   a0---b0---c0---d0   Flat plate with boundary layer
--
--             0---> x
mm = 1.0e-3 -- metres per mm
-- Leading edge of shock generator and inlet to the flow domain.
L1 = 10.0*mm; H1 = 37.36*mm
a0 = Vector3:new{x=-L1, y=0.0}
a1 = a0+Vector3:new{x=0.0,y=H1}
-- Angle of inviscid shock generator.
alpha = 3.09*math.pi/180.0
tan_alpha = math.tan(alpha)
-- Start of flat plate with boundary layer.
b0 = Vector3:new{x=0.0, y=0.0}
b1 = b0+Vector3:new{x=0.0,y=H1-L1*tan_alpha}
-- End of shock generator is only part way long the plate.
L3 = 67*mm
c0 = Vector3:new{x=L3, y=0.0}
c1 = c0+Vector3:new{x=0.0,y=H1-(L1+L3)*tan_alpha}
-- End of plate, and of the whole flow domain.
L2 = 90.0*mm
d0 = Vector3:new{x=L2, y=0.0}
d1 = d0+Vector3:new{x=0.0,y=H1}
-- Now, define the three patches.
patch0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
patch1 = CoonsPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}
patch2 = CoonsPatch:new{p00=c0, p10=d0, p11=d1, p01=c1}
--
-- Discretization of the flow domain.
--
-- We want to cluster the cells toward the surface of the flat plate.
-- where the boundary layer will be developing.
rcf = RobertsFunction:new{end0=true,end1=true,beta=1.1}
factor = 2 -- We'll scale discretization off this value
ni0 = math.floor(20*factor); nj0 = math.floor(80*factor)
grid0 = StructuredGrid:new{psurface=patch0, niv=ni0+1, njv=nj0+1,
			   cfList={east=rcf,west=rcf}}
grid1 = StructuredGrid:new{psurface=patch1, niv=7*ni0+1, njv=nj0+1,
			   cfList={east=rcf,west=rcf}}
grid2 = StructuredGrid:new{psurface=patch2, niv=2*ni0+1, njv=nj0+1,
			   cfList={east=rcf,west=rcf}}
--
-- Build the flow blocks and attach boundary conditions.
--
-- blk0 = FluidBlockArray{grid=grid0, initialState=inflow, nib=1, njb=2,
-- 		       bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
-- 			       north=WallBC_WithSlip:new{},
-- 			       south=WallBC_WithSlip:new{}}}
-- blk1 = FluidBlockArray{grid=grid1, initialState=inflow, nib=7, njb=2,
-- 		       bcList={south=WallBC_NoSlip_Adiabatic:new{},
-- 			       north=WallBC_WithSlip:new{}}}
-- blk2 = FluidBlockArray{grid=grid2, initialState=inflow, nib=2, njb=2,
-- 		       bcList={south=WallBC_NoSlip_Adiabatic:new{},
-- 			       north=WallBC_WithSlip:new{},
-- 			       east=OutFlowBC_FixedPT:new{p_outside=p_inf,
-- 							  T_outside=T_inf}}}

blk0 = FluidBlockArray{grid=grid0, initialState=inflow, nib=1, njb=2,
		       bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
			       north=WallBC_WithSlip:new{},
			       south=WallBC_WithSlip:new{}}}
blk1 = FluidBlockArray{grid=grid1, initialState=inflow, nib=7, njb=2,
		       bcList={south=WallBC_NoSlip_Adiabatic:new{},
			       north=WallBC_WithSlip:new{}}}
blk2 = FluidBlockArray{grid=grid2, initialState=inflow, nib=2, njb=2,
		       bcList={south=WallBC_NoSlip_Adiabatic:new{},
			       north=WallBC_WithSlip:new{},
			       east=OutFlowBC_Simple:new{}}}

identifyBlockConnections()

--------------------- MPI Setup ---------------------
mpiTasks = mpiDistributeBlocks(4, "load-balance")
config.spatial_deriv_from_many_points = false

config.gasdynamic_update_scheme = "euler"
config.flux_calculator = 'adaptive'
config.viscous = true
config.spatial_deriv_calc = 'divergence'
config.cfl_value = 0.5
config.max_time = 5.0*L2/u_inf -- time in flow lengths
config.max_step = 200000
config.dt_init = 1.0e-8
config.dt_plot = config.max_time/1
