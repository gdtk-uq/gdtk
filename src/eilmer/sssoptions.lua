-- Module for setting steady-state-solver options, required by prep.lua.
--
-- Authors: RJG and Kyle D.
--

-- Storage for steady-state solver settings
sssOptionsHidden = { -- hidden from user
   -- set defaults here

   temporal_integration_mode = 0,

   -- preconditioner settings
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   precondition_matrix_flux_calculator = "same as config.flux_calculator",
   -- some parameters to help reduce the cost of forming precondition matrix
   frozen_preconditioner_count = 1, -- how often the precondition matrix is updated
   start_preconditioning = 1, -- what iteration to start preconditioning on
   ilu_fill = 0, -- level of fill-in for ILU decomposition (use 0 for practical simulations)
   preconditioner_sigma = 1.0e-30,

   frozen_limiter_on_lhs = false,
   use_adaptive_preconditioner = false,
   use_physicality_check = false,
   physicality_check_theta = 0.2,
   use_line_search = false,
   inviscid_cfl = false,
   use_scaling = true, -- always good to scale the linear system
   use_complex_matvec_eval = false,-- use complex variable Frechet derivative

   -- general simulation settings
   number_pre_steps = 10,
   number_total_steps = 100,
   max_number_attempts = 3,
   stop_on_relative_global_residual = 1.0e-99,
   stop_on_absolute_global_residual = 1.0e-99,
   stop_on_mass_balance = -1.0,

   -- DPLU-SGS max number of subiterations (e.g. kmax)
   max_sub_iterations = 1,

   -- Restarted preconditioned FGMRES settings
   max_outer_iterations = 10, -- higher settings typically mean fast convergence of
                              -- an iteration BUT extra cost
   max_restarts = 10, -- same as above, but not as effective, keep between 5-10
   number_inner_iterations = 5,

   -- CFL ramp settings
   include_turb_quantities_in_residual = true,
   residual_based_cfl_scheduling = true,
   cfl_max = 1.0e8,
   cfl_min = 1.0e-02,
   cfl_schedule_length = 0,
   cfl_schedule_value_list = {},
   cfl_schedule_iter_list = {},

   -- Options for start-up phase (first order interpolation applied during this phase)
   number_start_up_steps = 5,
   LHSeval0 = 1,
   RHSeval0 = 1,
   cfl0 = 1.0,
   eta0 = 0.5, -- to what level of relative residual drop will we converge an iteration
   tau0 = 0.1, -- relative residual at which the CFL will be ramped
   sigma0 = 1.0e-8, -- perturbation (for complex-step choose a VERY small number)
   p0 = 0.75, -- parameter used in the routine to increase the CFL

   -- Options for inexact Newton phase
   -- this phase will use second order interpolation, if selected
   -- settings are same as above
   LHSeval1 = 2,
   RHSeval1 = 2,
   cfl1 = 10.0,
   tau1 = 0.1,
   sigma1 = 1.0e-8,
   p1 = 1.0,
   eta_strategy = "constant",
   eta1 = 0.5,
   eta1_max = 0.9,
   eta1_min = 0.01,
   eta_ratio_per_step = 0.9,
   gamma = 0.9,
   alpha = 2.0,
   limiter_freezing_residual_reduction = 1e-99,
   limiter_freezing_count = 50,
   -- Options related to writing out snapshots and diagnostics
   snapshots_count = 10,
   number_total_snapshots = 5,
   write_diagnostics_count = 20,
   write_loads_count = 20,

   __index = function (t, k)
      return sssOptionsHidden[k]
   end,
   __newindex = function (t, k, v)
      if sssOptionsHidden[k] == nil then
	 print(string.format("The field '%s' cannot be set in 'SteadyStateSolver' table.", k))
      else
	 sssOptionsHidden[k] = v
      end
   end,
   __call = function (_, t)
      for k, v in pairs(t) do
	 sssOptionsHidden.__newindex(t, k, v)
      end
   end
}

SteadyStateSolver = {}
setmetatable(SteadyStateSolver, sssOptionsHidden)

-- Storage for shape sensitivity calculator settings
sscOptionsHidden = { -- hidden from user
   -- set defaults here
   pseudotime = false,
   pseudotime_lhs_jacobian_order = 1,
   adjoint_precondition_matrix_order = 0,
   read_frozen_limiter_values_from_file = false,
   -- sensitivity parameters
   epsilon = 1.0e-30,
   -- GMRES parameters
   maxOuterIterations = 10,
   maxRestarts = 10,
   cfl0=1.0,
   eta = 0.1,
   stop_on_relative_global_residual = 1.0e-99,
   -- Bezier curve fit parameters
   tol_bezier_curve_fit = 1.0e-06,
   max_steps_bezier_curve_fit = 10000,
   -- user-defined file
   user_defined_objective_file = "dummy-obj-file.lua",

   __index = function (t, k)
      return sscOptionsHidden[k]
   end,
   __newindex = function (t, k, v)
      if sscOptionsHidden[k] == nil then
	 print(string.format("The field '%s' cannot be set in 'ShapeSensitivityCalculator' table.", k))
      else
	 sscOptionsHidden[k] = v
      end
   end,
   __call = function (_, t)
      for k, v in pairs(t) do
	 sscOptionsHidden.__newindex(t, k, v)
      end
   end
}

ShapeSensitivityCalculator = {}
setmetatable(ShapeSensitivityCalculator, sscOptionsHidden)

-- Storage for solid domain loose update settings
sdluOptionsHidden = { -- hidden from user
   -- set defaults here
   max_newton_iterations = 10,
   tolerance_newton_update = 1.0e-2,
   max_gmres_iterations = 10,
   tolerance_gmres_solve = 1.0e-3,
   perturbation_size = 1.0e-2,

   __index = function (t, k)
      return sdluOptionsHidden[k]
   end,
   __newindex = function (t, k, v)
      if sdluOptionsHidden[k] == nil then
	 print(string.format("The field '%s' cannot be set in 'SolidDomainLooseUpdate' table.", k))
      else
	 sdluOptionsHidden[k] = v
      end
   end,
   __call = function (_, t)
      for k, v in pairs(t) do
	 sdluOptionsHidden.__newindex(t, k, v)
      end
   end
}

SolidDomainLooseUpdate = {}
setmetatable(SolidDomainLooseUpdate, sdluOptionsHidden)


  
