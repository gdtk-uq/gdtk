-- cylFlare.lua
-- KD, 13-Aug-2016
-- Adapted from Eilmer3 version: PJ, 11-May-2013, 29-May-2013
-- Model of the CUBRC hollow cylinder with extended-flare experiment. 

config.title = "Double-cone, sharp nose."
print(config.title)
config.axisymmetric = true

-- Conditions to match those reported for CUBRC Run 14
p_inf = 18.55   -- Pa
u_inf = 2576.0  -- m/s
T_inf = 102.2   -- degree K
T_vib = 2711.0  -- degrees K (but we will ignore for ideal-gas)
nsp, nmodes = setGasModel('ideal-N2-gas-model.lua')
inflow  = FlowState:new{p=p_inf, velx=u_inf, vely=0.0, T=T_inf}
initial = FlowState:new{p=p_inf/5, velx=0.0, vely=0.0, T=T_inf}
T_wall = 295.8 -- degree K


mm = 1.0e-3 -- metres per mm
R = 6.35 -- radii
theta = 65.0 * math.pi/180.0

p = Vector3:new{x=(8.68+R)*mm, y=0.0}
e0 = Vector3:new{x=8.68*mm, y=0.0}
e1 = Vector3:new{x=(8.68+R*(1.0-math.cos(theta)))*mm, y=mm*R*math.sin(theta)}
--a0 = Vector3:new{x=0.0, y=0.0}
-- a1 = Vector3:new{x=0.0,y=5*mm} -- leading edge of domain
a0 = Vector3:new{x=3.308603587*mm,y=0.0} -- leading edge of domain
a1 = Vector3:new{x=5.888873328*mm, y=7.336956933*mm}
b0 = Vector3:new{x=92.08*mm,y=42.94*mm} -- junction between cones
b1 = Vector3:new{x=76*mm,y=61*mm} -- out in the free stream
c0 = Vector3:new{x=153.69*mm,y=130.925*mm} -- downstream-edge of second cone
c1 = Vector3:new{x=124*mm,y=181*mm} -- out in the free stream
d0 = Vector3:new{x=193.68*mm,y=130.925*mm} -- down-stream edge of domain
d1 = Vector3:new{x=193.68*mm,y=181*mm}

bezier_p1 = Vector3:new{x=3.308603587*mm,y=20.336956933*mm}
bezier_p2 = Vector3:new{x=70*mm,y=50*mm}

a0b1 = Bezier:new{points={a0, bezier_p1, bezier_p2, b1}}
a0b1 = ArcLengthParameterizedPath:new{underlying_path=a0b1}
--a0e0 = Line:new{p0=a0, p1=e0}
e0a0 = Line:new{p0=e0, p1=a0}
--b1b0 = Line:new{p0=b1, p1=b0}
b0b1 = Line:new{p0=b0, p1=b1}
nose_arc = Arc:new{p0=e0, p1=e1, centre=p}
line = Line:new{p0=e1, p1=b0}
e0b0 = Polyline:new{segments={nose_arc, line}}

rcfx0 = RobertsFunction:new{end0=true,end1=false,beta=1.26}
rcfy0 = RobertsFunction:new{end0=true,end1=false,beta=1.26}

rcfx = RobertsFunction:new{end0=true,end1=false,beta=1.2}
rcfy = RobertsFunction:new{end0=true,end1=false,beta=1.1}


ni0 = 120; nj0 = 40 -- We'll scale discretization off these values
factor = 2.0
ni0 = ni0*factor; nj0 = nj0*factor

--psurfCone1 = makePatch{north=b1b0, east=e0b0, south=a0e0, west=a0b1}
psurfCone1 = makePatch{north=a0b1, east=b0b1, south=e0b0, west=e0a0}
grdCone1 = StructuredGrid:new{psurface=psurfCone1,
			    cfList= {north=rcfx0,east=rcfy,south=rcfx,west=rcfy},
			    niv=ni0+1, njv=nj0+1}

grdCone2 = StructuredGrid:new{psurface=CoonsPatch:new{p00=b0,p10=c0,p11=c1,p01=b1},
			    cfList= {north=None,east=rcfy,south=None,west=rcfy},
			    niv=ni0+1, njv=nj0+1}   

grdCone3 = StructuredGrid:new{psurface=CoonsPatch:new{p00=c0,p10=d0,p11=d1,p01=c1},
			    cfList= {north=None,east=rcfy,south=None,west=rcfy},
			    niv=tonumber(ni0/2.0)+1, njv=nj0+1}   

-- [TODO] Can we remove this boundary condition fudge now? PJ 2017-04-30
-- create a special boundary condition for the no_slip_fixed_T wall that doesn't reference KOmegaWall
LaminarWallBC = WallBC_NoSlip_FixedT:new{Twall=T_wall}
table.remove(LaminarWallBC.preSpatialDerivAction, 5)

cone1 = FBArray:new{grid=grdCone1, nib=6, njb=2,
			initialState=inflow,
			bcList={north=InFlowBC_Supersonic:new{flowState=inflow},
				east=None,
				south=LaminarWallBC,
				west=None},
			label="cone1"}

cone2 = FBArray:new{grid=grdCone2, nib=6, njb=2,
			initialState=inflow,
			bcList={north=InFlowBC_Supersonic:new{flowState=inflow},
				east=None,
				south=LaminarWallBC,
				west=None},
			label="cone2"}  

cone3 = FBArray:new{grid=grdCone3, nib=2, njb=2,
			initialState=inflow,
			bcList={north=InFlowBC_Supersonic:new{flowState=inflow},
				east=OutFlowBC_FixedP:new{p_outside=p_inf/5.0},
				south=LaminarWallBC,
				west=None},
			label="out"}  

identifyBlockConnections()  

-- Do a little more setting of global data.
config.use_viscosity_from_cells=true
config.adjust_invalid_cell_data = true
config.max_invalid_cells = 10
config.viscous = true
config.flux_calculator = "adaptive"
config.spatial_deriv_calc = "divergence"
config.spatial_deriv_locn = "faces"
--config.include_ghost_cells_in_spatial_deriv_clouds = true -- for the usg
--config.unstructured_limiter = "barth" -- for the usg
config.gasdynamic_update_scheme = "classic-rk3"
config.cfl_value = 1.0
config.viscous_signal_factor = 0.4
-- The settling of the separation bubble will probably dominate.
config.max_time = 2.0*2.5e-3 -- long enough, looking at earlier simulations
config.max_step = 2000000
config.dt_init = 1.0e-12
config.dt_plot = 0.25e-3

-- sketch.xaxis(0.0, 0.250, 0.05, -0.010)
-- sketch.yaxis(0.0, 0.250, 0.05, -0.010)
-- sketch.window(0.0, 0.0, 0.250, 0.250, 0.05, 0.05, 0.25, 0.25)
