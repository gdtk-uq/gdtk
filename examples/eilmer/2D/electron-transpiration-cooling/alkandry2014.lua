-- alkandry2014.lua
-- Simulation of Alkandry 2014 geometry. Inviscid shock fitting simulation to
-- extract shock shape for viscous simulation.
--
-- Hicham Alkandry, Kyle Hanquist, and Iain D. Boyd. 
-- "Conceptual Analysis of Electron Transpiration Cooling 
-- for the Leading Edges of Hypersonic Vehicles", 
-- 11th AIAA/ASME Joint Thermophysics and Heat Transfer Conference, 
-- AIAA AVIATION Forum, (AIAA 2014-2674) 
-- https://doi.org/10.2514/6.2014-2674 
--
-- Zachary J. Denman & Will O. Landsberg
-- 2017-12-11 & 2018-09-25

-- ===========================================================================
-- Inputs
-- ===========================================================================
-------------------
--    General    --
-------------------
job_title = "Alkandry 2014 - hypersonic leading edge."
print(config.title)
config.dimensions = 2
config.axisymmetric = false

-------------------
-- Time Stepping --
-------------------
N = 20  -- Number of flow lengths
N_solutions = 20
flow_length = 0.02 / 6000 -- 2 radaii / flow velocity
config.max_time = N*flow_length
config.dt_plot = config.max_time/N_solutions
config.dt_init = 1.0e-12
config.cfl_value = 0.5
config.max_step = 100000000

config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "euler"
---------------------
-- Viscous Effects --
---------------------
config.viscous = true
config.viscous_signal_factor = 0.2

-----------------------
-- Loads calculation --
-----------------------
config.write_loads = true
config.dt_loads = config.max_time/N_solutions
config.boundary_groups_for_loads = "cylinder-wedge"
-- config.boundary_groups_for_loads = "stagnation-line"

-------------------
-- Thermoionic emmision delay
-------------------
-- Activate thermionic energy balance after flow has established
-- Initially runs as adiabatic wall.
-- Speeds up convergence
config.thermionic_emission_bc_time_delay = 2*flow_length 

-- ===========================================================================
-- Gas and Flow conditions
-- ===========================================================================
-- Gas model
nsp, nmodes, gmodel = setGasModel('thermally-perfect-air-11sp.lua')

T_inf = 238.0  -- K
rho_inf = 2.3e-4  -- kg/m^3
mass_fraction = {N2=0.767, O2=0.233}

-- Use the gas package to compute the freestream pressure so that when we
-- change gas models everything continues to work (not that we would see very
-- much change anyway).
Q = GasState:new{gmodel}
Q.T = T_inf;
Q.rho = rho_inf;
Q.massf = mass_fraction
gmodel:updateThermoFromRHOT(Q); gmodel:updateSoundSpeed(Q)
p_inf = Q.p  -- Pa
u_inf = 6000.0  -- m/s
M_inf = u_inf/Q.a

-- Create FlowState objects
initial = FlowState:new{p=p_inf/5, T=T_inf, velx=0.0, vely=0.0, massf=mass_fraction}
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0, massf=mass_fraction}

-- ===========================================================================
-- Geometry - looks complicated, but the grid is quite simple
-- ===========================================================================
-- Translate the imported shock fitted curve a bit to capture the shock
translation_distance = -0.003
R = 0.010  -- m, nose radius
L = 2*R + translation_distance  -- m, body length, Alkandry was 0.250 m, shortened due to compute
theta = 5  -- deg, half angle of wedge, Alkandry (2014) use alpha
theta_rad = math.rad(theta)  -- rad
-- Points
Ax = 0.0; Ay = 0.0  -- Stagnation point
Px = Ax+R; Py = Ay  -- Centre of blunt nose
Bx = Px+R*math.cos(math.pi/2+theta_rad); By = Py+R*math.sin(math.pi/2+theta_rad)
Cx = Ax+L; Cy = By+(L-Bx)*math.tan(theta_rad)
-- Create Vector3 objects
p = Vector3:new{x=Px, y=Py}
q = Vector3:new{x=Qx, y=Qy}
a = Vector3:new{x=Ax, y=Py}
b = Vector3:new{x=Bx, y=By}
c = Vector3:new{x=Cx, y=Cy}
-- Lines
ab = Arc:new{p0=a, p1=b, centre=p}  -- Nose tip
bc = Line:new{p0=b, p1=c}  -- Wedge
abc =  Polyline:new{segments={ab, bc}}  -- Cylinder-wedge
ef = Spline2:new{filename="shock-shape-t0049-extracted.dat"}
gh = TranslatedPath:new{original_path=ef, shift=Vector3:new{x=translation_distance}}
-- Parametric surfaces, grids, and blocks
tpsrf1 = ChannelPatch:new{north=gh, south=abc, ruled=true, pure2D=true}
tpsrf2 = ChannelPatch:new{north=gh, south=abc, ruled=true, pure2D=true}
-- NORTH AND SOUTH BOUNDARIES
hc = ReversedPath:new{underlying_path=tpsrf1:make_bridging_path(1.0)}
ga = ReversedPath:new{underlying_path=tpsrf2:make_bridging_path(0.0)}
-- Full surface with north and south boundaries (just for info, not needed)
-- i.e., this is the orientation that is used later.
-- tpsrf = CoonsPatch:new{north=hc, south=ga, east=abc, west=gh}

-- Now split the full parametric surface up into some nicer pieces for better
-- orthogonality. We do this by creating some SubRangedPath objects and then
-- creating some new patches from those.
-- We only need to split up abc and gh at reasonable points so that we get a
-- reasonable level of orthogonality. We then use ChannelPatch to get some
-- nice boundaries in between to new Paths.
abc_split = 0.265
gh_split = 0.26
abc0 = SubRangedPath:new{underlying_path=abc, t0=0.0, t1=abc_split}
abc1 = SubRangedPath:new{underlying_path=abc, t0=abc_split, t1=1.0}
gh0 = SubRangedPath:new{underlying_path=gh, t0=0.0, t1=gh_split}
gh1 = SubRangedPath:new{underlying_path=gh, t0=gh_split, t1=1.0}

-- Use ChannelPatch to generate some straight boundaries
tpsrf3 = ChannelPatch:new{north=abc0, south=gh0, ruled=true, pure2D=true}
tpsrf4 = ChannelPatch:new{north=abc1, south=gh1, ruled=true, pure2D=true}

-- Extract the boundary between the two ChannelPatch(es)
-- This is the WEST boundary which we use as a NORTH boundary
mid = tpsrf3:make_bridging_path(1.0)

-- Final Surfaces
psrf0 = CoonsPatch:new{north=mid, south=ga, east=abc0, west=gh0}
psrf1 = CoonsPatch:new{north=hc, south=mid, east=abc1, west=gh1}

-- Clustering
-- cl_wall = RobertsFunction:new{end0=false, end1=true, beta=1.3}
cl_wall = RobertsFunction:new{end0=false, end1=true, beta=1.1}
cl_stag = RobertsFunction:new{end0=true, end1=false, beta=1.1}
cflist0 = {north=cl_wall, south=cl_wall, east=cl_stag, west=cl_stag}
cflist1 = {north=cl_wall, south=cl_wall, east=none, west=none}

-- Grids
wall_normal_cells = 30; stag_cells = 40; trailing_edge_cells = 40;
-- wall_normal_cells = 20; stag_cells = 20; trailing_edge_cells = 20;
grd0 = StructuredGrid:new{psurface=psrf0, niv=wall_normal_cells, njv=stag_cells, cfList=cflist0}
grd1 = StructuredGrid:new{psurface=psrf1, niv=wall_normal_cells, njv=trailing_edge_cells, cfList=cflist1}

-- Blocks
blk0 = FluidBlockArray{grid=grd0,
                       initialState=initial,
                       bcList={south=WallBC_WithSlip:new{group="stagnation-line"},
                               east=WallBC_ThermionicEmission:new{emissivity=1.0, Ar=1.20e6, 
                               phi=2, ThermionicEmissionActive=1, Twall_iterations=200, 
                               Twall_subiterations=50, group="cylinder-wedge"},
                               -- east=WallBC_NoSlip_Adiabatic:new{group="cylinder-wedge"},
                               west=InFlowBC_Supersonic:new{flowState=inflow}},
                       nib=2, njb=2}
blk1 = FluidBlockArray{grid=grd1,
                      initialState=initial,
                      bcList={north=OutFlowBC_Simple:new{},
                              east=WallBC_ThermionicEmission:new{emissivity=1.0, Ar=1.20e6, 
                              phi=2, ThermionicEmissionActive=1, Twall_iterations=200, 
                              Twall_subiterations=50, group="cylinder-wedge"},
                               -- east=WallBC_NoSlip_Adiabatic:new{group="cylinder-wedge"},
                              west=InFlowBC_Supersonic:new{flowState=inflow}},
                      nib=2, njb=2}

identifyBlockConnections()


-- ===========================================================================
-- Sketching stuff
-- ===========================================================================
-- Analytical estimates
delta = R*0.386*math.exp(4.67/M_inf^2)

print("Estimated shock stand-off distance= ", delta)
Rc = R*1.386*math.exp(1.8/(M_inf - 1)^0.75)

-- Points
Ax = 0.0; Ay = 0.0  -- Stagnation point
Px = Ax+R; Py = Ay  -- Centre of blunt nose
Fx = Ax-delta; Fy = Ay  -- Upstream point of grid on stagnation line
Qx = Fx+Rc; Qy = 0.0  -- Centre of radius of curvature of shock shape
Bx = Px+R*math.cos(math.pi/2+theta_rad); By = Py+R*math.sin(math.pi/2+theta_rad)
Cx = Ax+L; Cy = By+(L-Bx)*math.tan(theta_rad)
Ex = Px+Rc*math.cos(math.pi/2+theta_rad); Ey = Py+Rc*math.sin(math.pi/2+theta_rad)
E1x = Qx+Rc*math.cos(135*math.pi/180); E1y = Qy+Rc*math.sin(135*math.pi/180)
F1x = Qx+Rc*math.cos(165*math.pi/180); F1y = Qy+Rc*math.sin(165*math.pi/180)
Dx = Cx; Dy = 3.0*Cy
dofile("sketch_domain.lua")
