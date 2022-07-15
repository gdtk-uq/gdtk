-- cyl.lua
-- Cylinder in Mach 10 flow with type-IV shock-shock interaction.
--
-- Modelling the flow discussed in the paper:
-- J.N. Moss, T. Pot, B. Chanetz and M. Lefebvre
-- DSMC Simulation of Shock/Shock Interactions: Emphasis on Type IV Interactions.
--
-- PJ 2022-06-29

config.title = "Cylinder with type IV shock-shock interaction."
print(config.title)
config.dimensions = 2
config.viscous = true
config.spatial_deriv_locn = "vertices"
config.spatial_deriv_calc = "divergence"
--
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
-- Compute the flowstates either side of the oblique shock.
M1 = 10.0
gs1 = GasState:new{gm}
gs1.p = 5.9 -- Pa
gs1.T = 52.5 -- K
gm:updateThermoFromPT(gs1)
gm:updateSoundSpeed(gs1)
V1 = M1 * gs1.a
fs1 = FlowState:new{p=gs1.p, T=gs1.T, velx=V1}
print("fs1=", fs1)
--
beta = math.rad(24.9) -- estimated from figure 6
gs2, theta, V2 = gasflow.theta_oblique(gs1, V1, beta)
print("gs2=", fs2)
print("theta=", theta, "radians")
print("V2=", V2)
fs2 = FlowState:new{p=gs2.p, T=gs2.T,
                    velx=V2*math.cos(theta), vely=V2*math.sin(theta)}
print("fs2=", fs2)
--
initial = FlowState:new{p=gs1.p/5, T=gs1.T, velx=V1/2}
--
-- Build the domain as two quadrants about the cylinder.
local billig_patch = require "billig_patch"
R = 0.008 -- m
nxc = 60; nyc = 60; factor = 1
-- quadant above x-axis
bp = billig_patch.make_patch{Minf=M1, R=R, xc=R, x_scale=2.4, y_scale=3.0, quadrant=2}
cfw = RobertsFunction:new{end0=true, end1=false, beta=1.1}
cfe = RobertsFunction:new{end0=true, end1=false, beta=1.2}
cf = RobertsFunction:new{end0=false, end1=true, beta=1.2}
grid_q2 = StructuredGrid:new{psurface=bp.patch, niv=nxc*factor+1, njv=nyc*factor+1,
                             cfList={west=cfw, east=cfe, north=cf, south=cf}}
-- quadrant below x-axis
bp = billig_patch.make_patch{Minf=M1, R=R, xc=R, x_scale=2.4, y_scale=3.0, quadrant=3}
cfw = RobertsFunction:new{end0=false, end1=true, beta=1.1}
cfe = RobertsFunction:new{end0=false, end1=true, beta=1.2}
grid_q3 = StructuredGrid:new{psurface=bp.patch, niv=nxc*factor+1, njv=nyc*factor+1,
                             cfList={west=cfw, east=cfe, north=cf, south=cf}}
-- Combine
grid_q3:joinGrid(grid_q2, "north")
-- Flow Blocks and boundary conditions.
inflow_bc = InOutFlowBC_DualState:new{flowState1=fs1, flowState2=fs2,
                                      p={x=0.008,y=0.00814},
                                      n={x=-math.sin(beta),y=math.cos(beta)}}
surface_bc = WallBC_NoSlip_FixedT:new{Twall=300.0, group='loads'}
outflow_bc = OutFlowBC_Simple:new{}
blks = FBArray:new{grid=grid_q3, initialState=initial, label="blk", nib=1, njb=6,
                   bcList={west=inflow_bc, north=outflow_bc, south=outflow_bc, east=surface_bc}}
--
config.max_time = 600.0e-6  -- seconds
config.max_step = 400000
config.dt_plot = 20.0e-6
config.dt_history = 1.0e-6
config.dt_init = 1.0e-10
config.write_loads = true
config.dt_loads = 10.0e-6
config.viscous_signal_factor = 0.2
