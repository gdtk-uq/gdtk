-- The Ghia-Ghia-Shin shear-driven cavity, with near-incompressible flow.
--
-- Author: Peter J. and Rowan G.
-- Date: 2016-11-06
--
-- Ported from Eilmer3 example prepared by Peter J. [2014-03-17 -- 2014-05-26]
--
-- Ported to Newton-Krylov accelerator by Kyle A. Damm [2020-05-14]
--
-- Reference:
-- Ghia, U., Ghia, K. and Shin, C. (1982)
-- High-Re solutions for incompressible flow using the Navier-Stokes equations
-- and a multigrid method.
-- Journal of Computational Physics, 48:p.387
--

config.title = "Shear-driven cavity flow Ghia et al. 1982, Re=100"
config.dimensions = 2

setGasModel('ideal-air-gas-model.lua')
Re = 100
mu = 1.847e-5 -- Pa.s
uLid = 1.0 -- m/s
L = 1.0 -- m
rho = Re * mu / (uLid * L)
T = 300.0
print("rho= ", rho, " kg/m**3")
p = rho * 287.1 * T -- Pa
print("p= ", p, " Pa")

-- Start with zero velocity everywhere
initialFlow = FlowState:new{p=p, T=T, velx=0.0, vely=0.0}

-- Set up geometry
p00 = Vector3:new{x=0.0, y=0.0}; p10 = Vector3:new{x=L, y=0.0}
p01 = Vector3:new{x=0.0, y=L};   p11 = Vector3:new{x=L, y=L}

ncells = 129

-- Set up grid and blocks
clustering = RobertsFunction:new{end0=false, end1=false, beta=1.1}
grid = StructuredGrid:new{psurface=CoonsPatch:new{p00=p00, p10=p10, p11=p11, p01=p01},
			  niv=ncells+1, njv=ncells+1,
			  cfList={north=clustering, east=clustering,
				  south=clustering, west=clustering}
}

blk = FluidBlockArray{grid=grid, fillCondition=initialFlow, nib=2, njb=2,
		 bcList={north=WallBC_TranslatingSurface_Adiabatic:new{v_trans={x=uLid, y=0.0, z=0.0}},
			 east=WallBC_NoSlip_Adiabatic:new{},
			 south=WallBC_NoSlip_Adiabatic:new{},
			 west=WallBC_NoSlip_Adiabatic:new{}
			}
}

identifyBlockConnections()

-- convert to unstructured blocks
for i=1,4 do
     SBlock2UBlock(fluidBlocks[i])
end

-- convective flux settings
config.apply_limiter = false
config.extrema_clipping = false
config.interpolation_order = 2
config.flux_calculator = "ausm_plus_up"
config.freeze_limiter_on_step = 4000000

-- viscous flux settings
config.viscous = true
config.spatial_deriv_locn = "cells"
config.spatial_deriv_calc = "least_squares"
config.diffuse_wall_bcs_on_init = false

config.print_count = 20

SteadyStateSolver{
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   ilu_fill = 0,
   frozen_preconditioner_count = 5,
   start_preconditioning = 1,
   
   use_scaling = true,
   use_complex_matvec_eval = true,
   
   number_pre_steps = 20,
   number_total_steps = 800,
   stop_on_relative_global_residual = 1.0e-12,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 40,
   max_restarts = 5,

   -- Settings for start-up phase
   number_start_up_steps = 100,
   sigma0 = 1.0e-30,
   cfl0 = 10.0,
   eta0 = 0.1,
   tau0 = 1.0,
   p0 = 1.0,

   -- Settings for inexact Newton phase
   sigma1 = 1.0e-30,
   cfl1 = 10.0,
   tau1 = 1.0,
   eta1 = 0.1,
   p1 = 1.0,
   eta_strategy = "constant",

   -- Settings control write-out
   snapshots_count = 20,
   number_total_snapshots = 5,
   write_diagnostics_count = 1,
   write_loads_count = 1000,
}
