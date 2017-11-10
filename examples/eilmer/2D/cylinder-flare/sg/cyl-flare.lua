-- cylFlare.lua
-- KD, 13-Aug-2016
-- Adapted from Eilmer3 version: PJ, 11-May-2013, 29-May-2013
-- Model of the CUBRC hollow cylinder with extended-flare experiment. 

config.title = "Hollow cylinder with extended flare."
print(config.title)
config.axisymmetric = true

-- Conditions to match those reported for CUBRC Run 14
p_inf = 31.88   -- Pa
u_inf = 2304.0  -- m/s
T_inf = 120.4   -- degree K
T_vib = 2467.0  -- degrees K (but we will ignore for ideal-gas)
nsp, nmodes = setGasModel('ideal-N2-gas-model.lua')
inflow  = FlowState:new{p=p_inf, velx=u_inf, vely=0.0, T=T_inf}
initial = FlowState:new{p=p_inf/5, velx=0.0, vely=0.0, T=T_inf}
T_wall = 295.2 -- degree K

mm = 1.0e-3    -- metres per mm
L1 = 101.7*mm  -- cylinder length
L2 = 220.0*mm  -- distance to end of flare
R1 = 32.5*mm
alpha = 30.0*math.pi/180.0 -- angle of flare
tan_alpha = math.tan(alpha)
a0 = Vector3:new{x=0.0, y=R1}; a1 = a0+Vector3:new{x=0.0, y=5*mm} -- leading edge of cylinder
b0 = Vector3:new{x=L1, y=R1}; b1 = b0+Vector3:new{x=-5*mm, y=20*mm} -- start flare
c0 = Vector3:new{x=L2, y=R1+tan_alpha*(L2-L1)}; c1 = c0+Vector3:new{x=0.0,y=25*mm} -- end flare

rcfx = RobertsFunction:new{end0=true,end1=false,beta=1.2}
rcfy = RobertsFunction:new{end0=true,end1=false,beta=1.1}
ni0 = 200; nj0 = 80 -- We'll scale discretization off these values
factor = 1.0
ni0 = ni0*factor; nj0 = nj0*factor*2.5

grdCyl = StructuredGrid:new{psurface=CoonsPatch:new{p00=a0,p10=b0,p11=b1,p01=a1},
			    cfList= {north=rcfx,east=rcfy,south=rcfx,west=rcfy},
			    niv=ni0+1, njv=nj0+1}

grdFlare = StructuredGrid:new{psurface=CoonsPatch:new{p00=b0,p10=c0,p11=c1,p01=b1},
			    cfList= {north=None,east=rcfy,south=None,west=rcfy},
			    niv=ni0+1, njv=nj0+1}   

cyl = FluidBlockArray{grid=grdCyl, nib=6, njb=2,
		      initialState=inflow,
		      bcList={north=InFlowBC_Supersonic:new{flowState=inflow},
			      east=None,
			      south=WallBC_NoSlip_FixedT:new{Twall=T_wall},
			      west=InFlowBC_Supersonic:new{flowState=inflow}},
		      label="cyl"}

Flare = FluidBlockArray{grid=grdFlare, nib=6, njb=2,
			initialState=initial,
			bcList={north=InFlowBC_Supersonic:new{flowState=inflow},
				east=OutFlowBC_FixedP:new{p_outside=p_inf/5.0},
				south=WallBC_NoSlip_FixedT:new{Twall=T_wall},
				west=None},
			label="flare"}  

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
