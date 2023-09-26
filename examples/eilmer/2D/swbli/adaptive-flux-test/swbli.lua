-- swbli.lua
-- Simulation of a shockwave/boundary layer interaction, based on T4
-- experiments by Eric Chang
-- https://doi.org/10.1017/jfm.2020.877
-- @author: NNG (13/06/22)

-- Notes: This example has been cooked up to not converge with ausmdv alone
config.dimensions = 2
config.viscous = true
config.flux_calculator = "adaptive_hanel_ausmdv"
config.turbulence_model = "spalart_allmaras_edwards"
config.thermo_interpolator = "rhop"

-- It also fails to converge with the discrete (PJ) shock detector, unless
-- freezing is applied.
config.shock_detector = "NNG"
config.strict_shock_detector = false -- required for non-discrete shock detection
-- Smaller absolute values means more sensitive detection. Unfortunately, sims are 
-- quite sensitive to this number. Somewhere between -0.3  and -0.01 was okay for
-- this case, but feel free to experiment.
config.compression_tolerance = -0.2 

-- Gas model and flow conditions Mach 7 Enthalpy
nsp, nmodes, gm = setGasModel('gm-air.lua')
p_inf = 2624.0 -- Pa
u_inf = 2221.0 -- m/s
T_inf = 250.0  -- K

-- Compute the gas state in physical units from its nondimensional description
gs = GasState:new{gm}
gs.T = T_inf
gs.p = p_inf
gs.massf = {air=1.0}
gm:updateThermoFromPT(gs)
gm:updateTransCoeffs(gs)

-- Use updated gas properties to estimate turbulence quantities
turb_lam_viscosity_ratio = 5.0 -- Fully turbulent, equation (8) from Allmaras (2012)
nu_inf = gs.mu/gs.rho
nuhat_inf = turb_lam_viscosity_ratio*nu_inf

inflow = FlowState:new{p=p_inf, velx=u_inf, T=T_inf, nuhat=nuhat_inf, massf={air=1.0}}

L = 320e-3
H = 80e-3
angle = math.rad(12.0)
s = 20e-3
G = 215e-3
dex = G*math.cos(angle)
dey = G*math.sin(angle)

a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=L, y=0.0}

c = Vector3:new{x=0.0, y=H}
d = c + s*Vector3:new{x=1.0, y=0.0}
e = d + Vector3:new{x=dex, y=-dey}
f = Vector3:new{x=L, y=e.y}

ab = Line:new{p0=a,p1=b}
ac = Line:new{p0=a,p1=c}
bf = Line:new{p0=b,p1=f}
cd = Line:new{p0=c,p1=d}
de = Line:new{p0=d,p1=e}
ef = Line:new{p0=e,p1=f}
cf = Polyline:new{segments={cd,de,ef}}

niv = 256; njv=64

patch = CoonsPatch:new{north=cf, east=bf, south=ab, west=ac}
cfx = RobertsFunction:new{end0=false,end1=true,beta=1.5}
cflist = {north=cfx, east=GeometricFunction:new{a=0.001, r=1.2, N=njv},
          south=cfx, west=GeometricFunction:new{a=0.002, r=1.2, N=njv}}
grid = StructuredGrid:new{psurface=patch, niv=niv, njv=njv, cfList=cflist}
blks = FBArray:new{grid=grid, nib=8, njb=2, fillCondition=inflow,
                   bcList={east=OutFlowBC_Simple:new{},
                           south=WallBC_NoSlip_FixedT:new{Twall=300.0, group="wall"},
                           west=InFlowBC_Supersonic:new{flowCondition=inflow}}}

identifyBlockConnections()

-- loads settings
config.boundary_groups_for_loads = "wall"
config.write_loads = true
config.with_local_time_stepping=true

SteadyStateSolver{
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   ilu_fill = 0,
   frozen_preconditioner_count = 20,
   
   use_scaling = true,
   use_complex_matvec_eval = true,
   use_physicality_check = true,
   
   number_total_steps = 1200,
   stop_on_relative_global_residual = 1.0e-9,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 40,
   max_restarts = 4,

   residual_based_cfl_scheduling = true,
   cfl_max = 1e4,

   -- Settings for start-up phase
   number_start_up_steps = 0,
   cfl0 = 1.0,
   eta0 = 0.01,
   sigma0 = 1.0e-50,

   -- Settings for inexact Newton phase
   cfl1 = 1.0,
   sigma1 = 1.0e-50,
   eta1 = 0.01,
   tau1 = 0.5,
   p1 = 0.7,

   -- Settings control write-out
   snapshots_count = 100,
   number_total_snapshots = 100,
   write_diagnostics_count = 1,
   write_loads_count = 100,
}
