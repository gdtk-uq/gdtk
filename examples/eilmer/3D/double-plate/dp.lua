-- 3D Flat Plate pair for testing corner cell issues.
-- It bears a remarkable resemblance, in both geometry and flow condition,
-- to a certain thermal compression experiment...
-- 2023-02-27: NNG
--
config.dimensions = 3
config.flux_calculator = "ausmdv"
config.viscous = true
config.turbulence_model = "spalart_allmaras_edwards"
config.with_local_time_stepping=true

config.unstructured_limiter = "svan_albada"
config.extrema_clipping=false                     
config.scale_species_after_reconstruction = false 
config.enforce_species_density_positivity = false 
config.smooth_limiter_coeff = 1e-3
config.suppress_reconstruction_at_boundaries = true -- IMPORTANT FOR THIS CASE

nsp, nmodes, gm = setGasModel('ideal-air.lua')
p_inf = 50.0e3  -- Pa
u_inf = 2600.0  -- m/s
T_inf = 900.0   -- K

-- Set up gas state and update thermodynamic transfer coefficients
gas_inf = GasState:new{gm}
gas_inf.p = p_inf; gas_inf.T = T_inf
gm:updateThermoFromPT(gas_inf); gm:updateSoundSpeed(gas_inf); gm:updateTransCoeffs(gas_inf)

-- Turbulence quantities estimate
turb_lam_viscosity_ratio = 5.0 -- Transitional starting ratio from LARC website
nu_inf = gas_inf.mu/gas_inf.rho
nuhat_inf = turb_lam_viscosity_ratio*nu_inf
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, nuhat=nuhat_inf}

-- Read in grid:
for blockfiles in io.popen("ls -1 ./su2grid | wc -l"):lines() do
   nblocks = tonumber(blockfiles)
end

grids = {}
for i=0,nblocks-1 do
   fileName = string.format("su2grid/block_%d_grid.su2", i)
   grids[i] = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1.0}
end

my_bcDict = {inflow=InFlowBC_Supersonic:new{flowCondition=inflow},
             wall=WallBC_NoSlip_FixedT:new{Twall=300.0},
             outflow=OutFlowBC_Simple:new{},
             METIS_INTERIOR=ExchangeBC_MappedCell:new{cell_mapping_from_file=true}}

blks = {}
for i=0,nblocks-1 do
   blks[i] = FluidBlock:new{grid=grids[i], initialState=inflow, bcDict=my_bcDict}
end

SteadyStateSolver{
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   frozen_preconditioner_count = 50,
   
   use_scaling = true,
   use_complex_matvec_eval = true,
   use_physicality_check = true,
   
   number_total_steps = 900,
   stop_on_relative_global_residual = 1.0e-8,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 30,
   max_restarts = 0,
   residual_based_cfl_scheduling = true,
   cfl_max = 1e3,

   -- Settings for start-up phase
   number_start_up_steps = 0,
   cfl0 = 5e-1,
   tau0 = 0.01,
   eta0 = 0.01,
   sigma0 = 1.0e-30,
   p0 = 0.5,

   -- Settings for inexact Newton phase
   cfl1 = 2e0,
   tau1 = 0.01,
   eta1 = 0.01,
   sigma1 = 1.0e-30,
   p1 = 0.7,
   eta_strategy = "constant",

   -- Settings control write-out
   snapshots_count = 100,
   number_total_snapshots = 2,
   write_diagnostics_count = 1,
   write_loads_count = 1000,
}
