-- fp.lua : Transitional Turbulent flow over a flat plate
-- Nick. N. Gibbons, based on simulations by Yu Chen
-- Conditions from:
--  "Transition of compressible high enthalpy boundary layer flow over a flat plate"
--  Y. He and R. G. Morgan
--  Aeronautical Journal, February 1994

config.title = "Mach 6.5 flow over a flat plate (SA-BCM)"
print(config.title)
config.dimensions = 2
config.turbulence_model = "spalart_allmaras_bcm"
config.freestream_turbulent_intensity = 0.017 -- Chosen to match experiments
config.viscous = true
config.flux_calculator = "ausmdv"

-- Gas model and flow conditions to match Table 1, the first entry
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
p_inf = 3.25e3  -- Pa
u_inf = 2100.0  -- m/s
T_inf = 254     -- K

-- Set up gas state and update thermodynamic transfer coefficients
gas_inf = GasState:new{gm}
gas_inf.p = p_inf; gas_inf.T = T_inf
gm:updateThermoFromPT(gas_inf); gm:updateSoundSpeed(gas_inf); gm:updateTransCoeffs(gas_inf)

-- Turbulence quantities estimate
turb_lam_viscosity_ratio = 0.015 -- Transitional starting ratio from LARC website
nu_inf = gas_inf.mu/gas_inf.rho
nuhat_inf = turb_lam_viscosity_ratio*nu_inf

-- Set up flow state
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, nuhat=nuhat_inf}

-- Geometry of the flow domain
L = 600.0e-3 -- metres
H = 0.20 * L
--
--         wall
--        c---------b
-- flow=> |         |
--        d         |
--          -\-     |
--    flow=>    -\- |
--        0         a ----> x
--
a = Vector3:new{x=L, y=0.0}; b = Vector3:new{x=L, y=H};
c = Vector3:new{x=0.0, y=H}; d = Vector3:new{x=0.0, y=3.0*H/4.0}
patch = CoonsPatch:new{p00=d, p10=a, p11=b, p01=c}

-- For testing
niv = 64+1; njv = 48+1;
cflist = {north=cfx, east=GeometricFunction:new{a=0.0005, r=1.2, N=njv, reverse=true},
          south=cfx, west=GeometricFunction:new{a=0.0010, r=1.2, N=njv, reverse=true}}

-- For actual converged simulations
--niv = 128+1; njv = 64+1;
--cflist = {north=cfx, east=GeometricFunction:new{a=0.0002, r=1.2, N=njv, reverse=true},
--          south=cfx, west=GeometricFunction:new{a=0.0004, r=1.2, N=njv, reverse=true}}
cfx = RobertsFunction:new{end0=true,end1=false,beta=1.05}
grd = StructuredGrid:new{psurface=patch, niv=niv, njv=njv, cfList=cflist}

blks = FBArray:new{grid=grd, nib=4, njb=2, fillCondition=inflow,
                   bcList={north=WallBC_NoSlip_FixedT:new{Twall=300.0, group="wall"},
                           east=OutFlowBC_Simple:new{},
                           south=InFlowBC_Supersonic:new{flowCondition=inflow},
                           west=InFlowBC_Supersonic:new{flowCondition=inflow}}}

identifyBlockConnections()

-- loads settings
config.boundary_groups_for_loads = "wall"
config.write_loads = true

SteadyStateSolver{
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   ilu_fill = 0,
   frozen_preconditioner_count = 50,
   start_preconditioning = 0,
   
   use_scaling = true,
   use_complex_matvec_eval = true,
   
   number_total_steps = 6000,
   stop_on_relative_global_residual = 1.0e-8,

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
   tau1 = 10.0,

   -- Settings control write-out
   snapshots_count = 100,
   number_total_snapshots = 100,
   write_diagnostics_count = 1,
   write_loads_count = 100,
}
