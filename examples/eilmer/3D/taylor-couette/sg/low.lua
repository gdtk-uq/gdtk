-- Author: Kan Qin
-- Date: 2016-06-05
--
-- Taylor-Couette flow
--
-- Set simulation parameters
  -- title
config.title = "Taylor Couette flow."
print(config.title)
  -- dimensions
config.dimensions = 3
  -- inviscid flux
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.interpolate_in_local_frame = true
  -- viscous flux
config.viscous = true
  -- printing
config.print_count = 20
config.dt_init = 1.0e-12
  -- steady state solver
SteadyStateSolver{
   use_preconditioning = false,
   -- sigma = 1.0e-6, -- presently it's computed internally
   number_pre_steps = 10,
   number_total_steps = 2000,
   -- Settings for FGMRES iterative solver
   max_outer_iterations = 60,
   -- number_inner_iterations = 5, -- not needed is preconditioning is false
   max_restarts = 9,
   -- Settings for start-up phase
   number_start_up_steps = 30,
   cfl0 = 1.0,
   eta0 = 0.5,
   tau0 = 0.1,
   sigma0 = 2.0e-6,
   -- Settings for inexact Newton phase
   cfl1 = 2.0,
   tau1 = 0.1,
   sigma1 = 1.0e-6,
   eta1 = 0.01,
   eta_strategy = "constant",
   -- Settings control write-out
   snapshots_count = 200,
   number_total_snapshots = 200,
   write_diagnostics_count = 1
}

-- Flow conditions, free stream
p_inf   = 100.0
r_omega = 2.0*math.pi*27600.0/60.0
T_1     = 348.0 
T_2     = 350.0
theta_A = 0.0*math.pi/180
theta_B = 5.0*math.pi/180 
theta   = theta_B - theta_A

-- Geometry
r_1     = 0.2125 
g_width = 0.0031 
r_2     = r_1 + g_width 
h_1     = 0.0 
h_2     = 10.0*g_width 

-- Grid dimensions
nx = 30
ny = 40
nz = 40

nsp, nmodes, gm = setGasModel('ideal-n2-gas-model.lua')

-- Set up inflow condition
inflow = FlowState:new{p=p_inf, T=(T_1+T_2)/2.0, velx=0.0, vely=0.0, velz=0.0}

-- set up computational domain
  -- points
a = Vector3:new{x=r_1*math.cos(theta_A), y=0.0*math.sin(theta_A), z=h_1}
b = Vector3:new{x=r_2*math.cos(theta_A), y=0.0*math.sin(theta_A), z=h_1}
c = Vector3:new{x=r_2*math.cos(theta_B), y=r_2*math.sin(theta_B), z=h_1}
d = Vector3:new{x=r_1*math.cos(theta_B), y=r_1*math.sin(theta_B), z=h_1}
e = Vector3:new{x=r_1*math.cos(theta_A), y=0.0*math.sin(theta_A), z=h_2}
f = Vector3:new{x=r_2*math.cos(theta_A), y=0.0*math.sin(theta_A), z=h_2}
g = Vector3:new{x=r_2*math.cos(theta_B), y=r_2*math.sin(theta_B), z=h_2}
h = Vector3:new{x=r_1*math.cos(theta_B), y=r_1*math.sin(theta_B), z=h_2}
centre0 = Vector3:new{x=0.0, y=0.0, z=h_1}
centre1 = Vector3:new{x=0.0, y=0.0, z=h_2}
  -- lines
ab = Line:new{p0=a, p1=b}
bc = Arc:new{p0=b, p1=c, centre=centre0}
dc = Line:new{p0=d, p1=c}
ad = Arc:new{p0=a, p1=d, centre=centre0}
ef = Line:new{p0=e, p1=f}
fg = Arc:new{p0=f, p1=g, centre=centre1}
hg = Line:new{p0=h, p1=g}
eh = Arc:new{p0=e, p1=h, centre=centre1}
ae = Line:new{p0=a, p1=e}
bf = Line:new{p0=b, p1=f}
cg = Line:new{p0=c, p1=g}
dh = Line:new{p0=d, p1=h}
  -- Create a surface from this geometry in the xy-plane
surf_n = makePatch{north=hg, east=cg, south=dc, west=dh}
surf_e = makePatch{north=fg, east=cg, south=bc, west=bf}
surf_s = makePatch{north=ef, east=bf, south=ab, west=ae}
surf_w = makePatch{north=eh, east=dh, south=ad, west=ae}
surf_b = makePatch{north=dc, east=bc, south=ab, west=ad}
surf_t = makePatch{north=hg, east=fg, south=ef, west=eh}
  -- Now set-up volume
clusterx = RobertsFunction:new{end0=true, end1=true, beta=1.1}
clustery = RobertsFunction:new{end0=true, end1=true, beta=1.1}
clusterz = RobertsFunction:new{end0=true, end1=true, beta=1.1}
cflist = {edge01=clusterx, edge12=clustery, edge32=clusterx, edge03=clustery,
	  edge45=clusterx, edge56=clustery, edge76=clusterx, edge47=clustery,
	  edge04=clusterz, edge15=clusterz, edge26=clusterz, edge37=clusterz}
vol = TFIVolume:new{north=surf_n, east=surf_e, south=surf_s, west=surf_w, top=surf_t, bottom=surf_b}
  -- Set up grid
grid = StructuredGrid:new{pvolume=vol, cfList=cflist, niv=nx+1, njv=ny+1, nkv=nz+1}
  -- Set up block
bsouth = ExchangeBC_FullFace:new{otherBlock=3, otherFace=north,
                                     orientation=0, reorient_vector_quantities=true,
                                     Rmatrix={9.96194698091745545199e-01, 8.71557427476581658699e-02, 0.00000000000000000000e+00, -8.71557427476581658699e-02, 9.96194698091745545199e-01, 0.00000000000000000000e+00, 0.00000000000000000000e+00, 0.00000000000000000000e+00, 1.00000000000000000000e+00}}

bnorth = ExchangeBC_FullFace:new{otherBlock=0, otherFace=south,
                                     orientation=0, reorient_vector_quantities=true,
                                     Rmatrix={9.96194698091745545199e-01, -8.71557427476581658699e-02, 0.00000000000000000000e+00, 8.71557427476581658699e-02, 9.96194698091745545199e-01, 0.00000000000000000000e+00, 0.00000000000000000000e+00, 0.00000000000000000000e+00, 1.00000000000000000000e+00}}
beast   = WallBC_NoSlip_FixedT:new{Twall=T_2}
bwest   = WallBC_RotatingSurface_FixedT:new{Twall=T_1, r_omega={x=0.0, y=0.0, z=r_omega}, centre={x=0.0, y=0.0, z=0.0}}
btop    = WallBC_WithSlip:new{}
bbottom = WallBC_WithSlip:new{}
bclist = {north=bnorth, east=beast, south=bsouth, west=bwest, top=btop, bottom=bbottom}
blk = FluidBlockArray{grid=grid, nib=1, njb=4, nkb=1, bcList=bclist, fillCondition=inflow, label="block-0"}

identifyBlockConnections()
