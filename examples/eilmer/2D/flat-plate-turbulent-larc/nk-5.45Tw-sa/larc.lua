-- larc.lua : Turbulent flow over a flat plate
-- 2D Zero Pressure Gradient High Mach Number Flat Plate Validation Case
-- https://turbmodels.larc.nasa.gov/ZPGflatplateSS_val_sa.html
-- @author: Nick N. Gibbons (n.gibbons@uq.edu.au)
config.title = "LARC Mach 5.0 flow over a flat plate (spalart-allmaras)"
print(config.title)
config.dimensions = 2
config.turbulence_model = "spalart_allmaras"
config.viscous = true
config.flux_calculator = "ausmdv"
config.interpolation_order = 2

-- Gas model and flow conditions to match ZPG case M=5, Tw/Tinf=5.450
-- Shout out to this incredibly inconvenient problem description...
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
Re_unit = 15e6
T_inf = 540.6*5/9 -- (The 540.6 is degrees Rankine, convert to Kelvin)
M_inf = 5.0
Tw_on_Tinf = 5.450
Tw = Tw_on_Tinf*T_inf

-- Compute the gas state in physical units from its nondimensional description
gas_inf = GasState:new{gm}
gas_inf.T = T_inf
gm:updateSoundSpeed(gas_inf)
u_inf = M_inf*gas_inf.a
gm:updateTransCoeffs(gas_inf)
gas_inf.rho = Re_unit*gas_inf.mu/u_inf
gm:updateThermoFromRHOT(gas_inf)
p_inf = gas_inf.p

-- Use updated gas properties to estimate turbulence quantities
turb_lam_viscosity_ratio = 5.0 -- Fully turbulent, equation (8) from Allmaras (2012)
nu_inf = gas_inf.mu/gas_inf.rho
nuhat_inf = turb_lam_viscosity_ratio*nu_inf

inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, nuhat=nuhat_inf}
--inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf}
print("Inflow Check\n", inflow)
print("gas_inf.rho: ", gas_inf.rho, "gas_inf.mu: ", gas_inf.mu)


-- Set up grid
grd0 = StructuredGrid:new{filename="../grid_35x25/larc.b0000.txt", fmt="text"}
grd1 = StructuredGrid:new{filename="../grid_35x25/larc.b0001.txt", fmt="text"}

--  *--- SupIn ---* ------ SupIn --------*
--  |             |                      |
--  |             |                      |
--SupIn->  0      |           1          Out ->
--  |             |                      |
--  |             |                      |
--  *--Symmetry---*------FixedTWall------*

blk0 = FBArray:new{grid=grd0, initialState=inflow, nib=1, njb=2, 
            bcList={north=InFlowBC_Supersonic:new{flowCondition=inflow},
                    west=InFlowBC_Supersonic:new{flowCondition=inflow}}
}

blk1 = FBArray:new{grid=grd1, initialState=inflow, nib=2, njb=2,
          bcList={south=WallBC_NoSlip_FixedT0:new{Twall=Tw,wall_function=false,group="wall"},
--            bcList={south=WallBC_NoSlip_Adiabatic:new{wall_function=false,group="wall"},
                    north=InFlowBC_Supersonic:new{flowCondition=inflow},
                    east=OutFlowBC_Simple:new{}}
}

identifyBlockConnections()

for i=1,6 do
     SBlock2UBlock(fluidBlocks[i])
end

-- loads settings
config.boundary_groups_for_loads = "wall"
config.write_loads = true
config.freeze_limiter_on_step = 1400

-- Steady State Settings, to match Mabey example case
config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"
config.diffuse_wall_bcs_on_init = true
config.number_init_passes = 25

SteadyStateSolver{
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   ilu_fill = 0,
   frozen_preconditioner_count = 25,
   start_preconditioning = 1,
   
   use_scaling = true,
   use_complex_matvec_eval = true,
   
   number_pre_steps = 10,
   number_total_steps = 2000,
   stop_on_relative_global_residual = 1.0e-13,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 30,
   max_restarts = 10,

   -- Settings for start-up phase
   number_start_up_steps = 200,
   cfl0 = 0.5,
   eta0 = 0.1,
   tau0 = 0.5,
   sigma0 = 1.0e-30,
   p0 = 0.5,

   -- Settings for inexact Newton phase
   cfl1 = 10.0,
   tau1 = 1.0,
   sigma1 = 1.0e-30,
   eta1 = 0.01,
   p1 = 1.0,
   eta1_min = 0.01,
   eta_ratio_per_step = 0.99,
   eta_strategy = "geometric",

   -- Settings control write-out
   snapshots_count = 2000,
   number_total_snapshots = 100,
   write_diagnostics_count = 1,
   write_loads_count = 10000,
}
