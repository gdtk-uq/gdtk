-- configoptions.lua
-- Module for setting transient solver options, required by prep.lua.
--
-- Authors: PAJ, RJG, Kyle D, Nick G, Daryl B and others.
--

local configoptions = {}

configOptionsHidden = { -- hidden from user
   -- set defaults here
   base_file_name = "job",
   grid_format = "gziptext",
   flow_format = "gziptext",
   new_flow_format = false,
   title = "Eilmer4 simulation",
   --
   gas_model_file = "gas-model.lua",
   sticky_electrons = false,
   include_quality = false,
   --
   udf_supervisor_file = "",
   user_pad_length = 0,
   nFluidBlocks = 0,
   dimensions = 2,
   true_centroids = false,
   axisymmetric = false,
   gravity = {x=0.0, y=0.0, z=0.0},
   --
   MHD = false,
   MHD_static_field = false,
   MHD_resistive = false,
   divergence_cleaning = false,
   c_h = 0.0,
   divB_damping_length = 1.0,
   electric_field_count = 1000000000,
   solve_electric_field = false,
   field_conductivity_model="none",
   --
   strang_splitting = "full_T_full_R",
   gasdynamic_update_scheme = "predictor_corrector",
   eval_udf_source_terms_at_each_stage = false,
   residual_smoothing = false,
   residual_smoothing_weight = 0.2,
   residual_smoothing_type = "explicit",
   residual_smoothing_iterations = 2,
   with_local_time_stepping = false,
   local_time_stepping_limit_factor = 10000,
   with_super_time_stepping_flexible_stages = false,
   max_attempts_for_step = 3,
   perturbation_for_real_differences = 1.0e-6,
   --
   solid_domain_cfl = 0.85,
   coupling_with_solid_domains = "tight",
   solid_has_isotropic_properties = true,
   solid_has_homogeneous_properties = true,
   solid_domain_augmented_deriv_avg = true,
   fluid_solid_bc_use_heat_transfer_coeff = false,
   --
   apply_bcs_in_parallel = true,
   --
   -- See struct FlowStateLimits in globalconfig.d.
   flowstate_limits_max_velocity = 30000.0, -- m/s
   flowstate_limits_max_tke = 0.01*1.0e38, -- guess for huge
   flowstate_limits_min_tke = 0.0,
   flowstate_limits_max_temp = 50000.0, -- Kelvin
   flowstate_limits_min_temp = 0.0, -- Kelvin
   flowstate_limits_min_pressure = 0.1, -- Pascals
   --
   user_specified_velocities_are_in_non_rotating_frame = true,
   --
   ignore_low_T_thermo_update_failure = true,
   suggested_low_T_value = 200.0,
   adjust_invalid_cell_data = false,
   report_invalid_cells = true,
   max_invalid_cells = 0,
   --
   n_ghost_cell_layers = 2,
   high_order_flux_calculator = false,
   flux_calculator = "adaptive_hanel_ausmdv",
   interpolation_order = 2,
   interpolation_delay = 0.0,
   allow_interpolation_for_sensitivity_matrix = false,
   suppress_radial_reconstruction_at_xaxis = false,
   suppress_reconstruction_at_shocks = false,
   suppress_reconstruction_at_boundaries = false,
   --
   thermo_interpolator = "rhou",
   allow_reconstruction_for_species = true,
   allow_reconstruction_for_energy_modes = true,
   allow_reconstruction_for_turbulent_variables = true,
   apply_limiter = true,
   epsilon_van_albada = 1e-12,
   extrema_clipping = true,
   apply_heuristic_pressure_based_limiting = false,
   interpolate_in_local_frame = true,
   apply_entropy_fix = true,
   enforce_species_density_positivity = false,
   scale_species_after_reconstruction = true,
   unstructured_limiter = "venkat",
   freeze_limiter_on_step = 1000000000,
   use_extended_stencil = false,
   smooth_limiter_coeff = 0.3,
   nsteps_of_chemistry_ramp = -1,
   shear_tolerance = 0.20,
   M_inf = 0.01,
   --
   shock_detector = "PJ",
   compression_tolerance = -0.30,
   do_shock_detect = false,
   damped_outflow = false,
   strict_shock_detector = true,
   shock_detector_smoothing = 0,
   frozen_shock_detector = false,
   shock_detector_freeze_step = 1000000000,
   --
   artificial_compressibility = false,
   ac_alpha = 0.1,
   --
   radiation = false,
   --
   grid_motion = "none",
   write_vertex_velocities = false,
   udf_grid_motion_file = "dummy-grid-motion-file.txt",
   --
   shock_fitting_delay = 0.0,
   shock_fitting_allow_flow_reconstruction = true,
   shock_fitting_scale_factor = 0.5,
   shock_fitting_filter_velocity_scale = 0.0,
   shock_fitting_assume_symmetry_at_first_point = false,
   --
   viscous = false,
   use_viscosity_from_cells = false,
   spatial_deriv_calc = "least_squares",
   spatial_deriv_locn = "cells",
   include_ghost_cells_in_spatial_deriv_clouds = true,
   upwind_vertex_gradients = true,
   save_convective_gradients = false,
   save_viscous_gradients = false,
   save_limiter_values = false,
   save_residual_values = false,
   save_timestep_values = false,
   nic_write = 1,
   njc_write = 1,
   nkc_write = 1,
   viscous_factor_increment = 0.01,
   viscous_delay = 0.0,
   shear_stress_relative_limit = 1.0,
   apply_shear_stress_relative_limit = false,
   viscous_signal_factor = 1.0,
   turbulent_signal_factor = 1.0,
   mass_diffusion_model = "none",
   diffusion_coefficient_type = "none",
   lewis_number = 1.0,
   --
   turbulence_model = "none",
   turbulence_prandtl_number = 0.89,
   turbulence_schmidt_number = 0.75,
   max_mu_t_factor = 3000.0,
   transient_mu_t_factor = 1.0,
   freestream_turbulent_intensity = 1.0,
   tci_model = "none",
   --
   chemistry_update = "split", -- or "integral"
   reacting = false,
   reactions_file = "chemistry.lua",
   reaction_time_delay = 0.0,
   reaction_fraction_schedule = {},
   T_frozen = 300.0,
   T_frozen_energy = 300.0,
   ignition_time_start = 0.0,
   ignition_time_stop = 0.0,
   energy_exchange_file = "energy-exchange.lua",
   radiation_energy_dump_allowed = false,
   radiation_energy_dump_temperature_limit = 30000.0,
   --
   max_step = 100,
   write_loads_at_step = -1,
   write_flow_solution_at_step = -1,
   snapshot_count = 1000000000,
   number_total_snapshots = 0,
   halt_now = 0,
   print_count = 20,
   control_count = 10,
   verbosity_level = 1,
   start_time = 0.0,
   max_time = 1.0e-3,
   --
   dt_init = 1.0e-3,
   dt_max = 1.0e-3,
   cfl_value = 0.5,
   -- If the user does not set a schedule of cfl values,
   -- we drop back to using the single cfl_value to construct the schedule.
   cfl_schedule = {},
   cfl_scale_factor = 1.0,
   stringent_cfl = false,
   cfl_count = 10,
   fixed_time_step = false,
   dt_plot = 1.0e-3,
   dt_history = 1.0e-3,
   --
   dt_loads = 1.0e-3,
   boundary_groups_for_loads = "loads",
   write_loads = false,
   compute_run_time_loads = false,
   run_time_loads_count = 100,
   --
   diffuse_wall_bcs_on_init = false,
   number_init_passes = 30,
   wall_temperature_on_init = -1.0,
   --
   block_marching = false,
   nib = 1,
   njb = 1,
   nkb = 1,
   propagate_inflow_data = false,
   save_intermediate_results = false,
   --
   udf_source_terms_file = "dummy-source-terms.txt",
   udf_source_terms = false,
   udf_solid_source_terms_file = "dummy-solid-source-terms.txt",
   udf_solid_source_terms = false,
   --
   do_temporal_DFT = false,
   DFT_n_modes = 5,
   DFT_step_interval = 10,
   --
   do_flow_average = false,
   --
   __index = function (t, k)
      return configOptionsHidden[k]
   end,
   __newindex = function (t, k, v)
      if configOptionsHidden[k] == nil then
	 print(string.format("The field '%s' cannot be set in 'config' table.", k))
      else
	 configOptionsHidden[k] = v
      end
   end,
   __call = function (_, t)
      for k, v in pairs(t) do
	 configOptionsHidden.__newindex(t, k, v)
      end
   end
} -- end table configOptionsHidden

configoptions.config = {}
setmetatable(configoptions.config, configOptionsHidden)

return configoptions
