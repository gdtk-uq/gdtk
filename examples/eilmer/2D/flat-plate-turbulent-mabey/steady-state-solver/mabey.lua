-- mabey.lua : Turbulent flow over a flat plate
-- Dimir Y.X. Pot, Wilson Y.K. Chan, 2018-03-14
-- Ported from Eilmer3
--  Mabey test case (AGARDograph 223 - Test series 7402)
--  (Referenced from Fernholz & Finley (1977),
--  AGARDograph No. 223, "A critical compilation of
--  compressible turbulent boundary layer data.")
--
config.title = "Mabey Mach 4.5 flow over a flat plate (k-omega)"
print(config.title)
config.dimensions = 2
config.turbulence_model = "k_omega"

-- Gas model and flow conditions to match Mabey's data set 74021802
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
p_inf = 3.16e3  -- Pa
u_inf = 712.9   -- m/s
T_inf = 62.16   -- K

-- Set up gas state and update thermodynamic transfer coefficients
gas_inf = GasState:new{gm}
gas_inf.p = p_inf; gas_inf.T = T_inf
gm:updateThermoFromPT(gas_inf); gm:updateSoundSpeed(gas_inf); gm:updateTransCoeffs(gas_inf)

-- Turbulence quantities estimate
tke_inf = 1.5 * (u_inf * 0.01)^2
mu_t_inf = 1.0 * gas_inf.mu
omega_inf = gas_inf.rho * tke_inf / mu_t_inf

-- Set up flow state
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, tke=tke_inf, omega=omega_inf}
print("Inflow Check\n", inflow)

-- Geometry of the flow domain
L = 1.40 -- metres
H = 0.40 * L
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
cfx = RobertsFunction:new{end0=true,end1=false,beta=1.05}
cflist = {north=cfx, east=RobertsFunction:new{end0=false,end1=true,beta=1.0014},
	  south=cfx, west=RobertsFunction:new{end0=false,end1=true,beta=1.0074}}
grd = StructuredGrid:new{psurface=patch, niv=129/3, njv=97/3, cfList=cflist}

blks = FBArray:new{grid=grd, nib=2, njb=2, fillCondition=inflow,
                   bcList={north=WallBC_NoSlip_Adiabatic0:new{wall_function=false, group="wall"},
                           east=OutFlowBC_Simple:new{},
                           south=InFlowBC_Supersonic:new{flowCondition=inflow},
                           west=InFlowBC_Supersonic:new{flowCondition=inflow}}}

identifyBlockConnections()

-- convert to unstructured blocks
for i=1,4 do
     SBlock2UBlock(fluidBlocks[i])
end

-- loads settings
config.boundary_groups_for_loads = "wall"
config.write_loads = true

-- convective flux settings
config.interpolation_order = 2
config.flux_calculator = "ausmdv"
config.freeze_limiter_on_step = 4000

-- viscous flux settings
config.viscous = true
config.spatial_deriv_locn = "cells"
config.spatial_deriv_calc = "least_squares"
config.diffuse_wall_bcs_on_init = false
config.number_init_passes = 25

SteadyStateSolver{
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   ilu_fill = 0,
   frozen_preconditioner_count = 25,
   start_preconditioning = 0,
   
   use_scaling = true,
   use_complex_matvec_eval = true,
   
   number_pre_steps = 10,
   number_total_steps = 13000,
   stop_on_relative_global_residual = 1.0e-12,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 50,
   max_restarts = 0,

   residual_based_cfl_scheduling = false,
   cfl_max = 1e6,
   cfl_schedule_length = 5,
   cfl_schedule_value_list = {1,  1e1, 1e2, 1e3, 1e4},
   cfl_schedule_iter_list =  {1,  20,  60,  85,  110},

   -- Settings for start-up phase
   number_start_up_steps = 0,
   cfl0 = 1,
   eta0 = 0.5,
   sigma0 = 1.0e-50,

   -- Settings for inexact Newton phase
   cfl1 = 1,
   sigma1 = 1.0e-50,
   eta1 = 0.5,

   -- Settings control write-out
   snapshots_count = 50,
   number_total_snapshots = 100,
   write_diagnostics_count = 1,
   write_loads_count = 10000,
}
