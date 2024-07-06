-- Module for writing files for preparing a simulation, required by prep.lua.
--
-- Authors: PJ, RJG, Kyle D. and Nick G.
--

local output = {}

function output.write_control_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"dt_init": %.18e,\n', config.dt_init))
   f:write(string.format('"dt_max": %.18e,\n', config.dt_max))
   f:write(string.format('"cfl_scale_factor": %.18e,\n', config.cfl_scale_factor))
   f:write(string.format('"stringent_cfl": %s,\n', tostring(config.stringent_cfl)))
   f:write(string.format('"viscous_signal_factor": %.18e,\n', config.viscous_signal_factor))
   f:write(string.format('"turbulent_signal_factor": %.18e,\n', config.turbulent_signal_factor))
   f:write(string.format('"residual_smoothing_weight": %.18e,\n', config.residual_smoothing_weight))
   f:write(string.format('"residual_smoothing_iterations": %d,\n', config.residual_smoothing_iterations))
   f:write(string.format('"residual_smoothing_type": "%s",\n',
			 config.residual_smoothing_type))
   f:write(string.format('"fixed_time_step": %s,\n', tostring(config.fixed_time_step)))
   f:write(string.format('"print_count": %d,\n', config.print_count))
   f:write(string.format('"cfl_count": %d,\n', config.cfl_count))
   f:write(string.format('"max_time": %.18e,\n', config.max_time))
   f:write(string.format('"max_step": %d,\n', config.max_step))
   f:write(string.format('"dt_plot": %.18e,\n', config.dt_plot))
   f:write(string.format('"dt_history": %.18e,\n', config.dt_history))
   f:write(string.format('"dt_loads": %.18e,\n', config.dt_loads))
   f:write(string.format('"write_loads_at_step": %d,\n', config.write_loads_at_step))
   f:write(string.format('"write_flow_solution_at_step": %d,\n', config.write_flow_solution_at_step))
   f:write(string.format('"snapshot_count": %d,\n', config.snapshot_count))
   f:write(string.format('"number_total_snapshots": %d,\n', config.number_total_snapshots))
   f:write(string.format('"halt_now": %d\n', config.halt_now))
   -- Note, also, no comma on last entry in JSON object. (^^^: Look up one line and check!)
   f:write("}\n")
   f:close()
end

function output.write_config_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"solver_mode": "%s",\n', config.solver_mode))
   f:write(string.format('"start_time": %.18e,\n', config.start_time))
   f:write(string.format('"grid_format": "%s",\n', config.grid_format))
   f:write(string.format('"field_format": "%s",\n', config.field_format))
   f:write(string.format('"gas_model_file": "%s",\n', config.gas_model_file))
   f:write(string.format('"udf_supervisor_file": "%s",\n', tostring(config.udf_supervisor_file)))
   if type(user_pad_data) == 'table' then
      if config.user_pad_length < #user_pad_data then
         config.user_pad_length = #user_pad_data
      end
   end
   f:write(string.format('"user_pad_length": %d,\n', config.user_pad_length))
   f:write('"user_pad_data": [')
   if type(user_pad_data) == 'table' then
      for i,e in ipairs(user_pad_data) do
         f:write(string.format('%.18e', e))
         if i < #user_pad_data then f:write(', ') end
      end
   end
   f:write('],\n')
   f:write(string.format('"sticky_electrons": %s,\n', tostring(config.sticky_electrons)))
   f:write(string.format('"include_quality": %s,\n', tostring(config.include_quality)))
   f:write(string.format('"dimensions": %d,\n', config.dimensions))
   f:write(string.format('"true_centroids": %s,\n', tostring(config.true_centroids)))
   f:write(string.format('"axisymmetric": %s,\n', tostring(config.axisymmetric)))
   config.gravity.x = config.gravity.x or 0.0
   config.gravity.y = config.gravity.y or 0.0
   config.gravity.z = config.gravity.z or 0.0
   f:write(string.format('"gravity": [%.18e, %.18e, %.18e],\n', config.gravity.x, config.gravity.y, config.gravity.z))
   f:write(string.format('"strang_splitting": "%s",\n', config.strang_splitting))
   f:write(string.format('"gasdynamic_update_scheme": "%s",\n', config.gasdynamic_update_scheme))
   f:write(string.format('"residual_smoothing": %s,\n', tostring(config.residual_smoothing)))
   f:write(string.format('"with_local_time_stepping": %s,\n', tostring(config.with_local_time_stepping)))
   f:write(string.format('"local_time_stepping_limit_factor": %d,\n', tostring(config.local_time_stepping_limit_factor)))
   f:write(string.format('"with_super_time_stepping": %s,\n', tostring(config.with_super_time_stepping)))
   f:write(string.format('"with_super_time_stepping_flexible_stages": %s,\n', tostring(config.with_super_time_stepping_flexible_stages)))
   f:write(string.format('"max_attempts_for_step": %d,\n', config.max_attempts_for_step))
   f:write(string.format('"ignore_low_T_thermo_update_failure": %s,\n', tostring(config.ignore_low_T_thermo_update_failure)))
   f:write(string.format('"suggested_low_T_value": %.18e,\n', config.suggested_low_T_value))
   f:write(string.format('"perturbation_for_real_differences": %.18e,\n', config.perturbation_for_real_differences))
   if #config.cfl_schedule > 0 then
      -- The table already have some values.
      -- We will presume that they are valid entries, however,
      -- if enough people get it wrong, we'll put some check here.
   else
      -- Fall back to making up a schedule from cfl_value.
      config.cfl_schedule = {{0.0, config.cfl_value},}
   end
   local cfl_schedule_length = #config.cfl_schedule
   f:write(string.format('"cfl_schedule_length": %d,\n', cfl_schedule_length))
   local cfl_values = {}
   local cfl_times = {}
   for i,cfl_pair in ipairs(config.cfl_schedule) do
      cfl_times[#cfl_times+1] = cfl_pair[1]
      cfl_values[#cfl_values+1] = cfl_pair[2]
   end
   f:write('"cfl_schedule_values": [')
   for i,e in ipairs(cfl_values) do
      f:write(string.format('%.18e', e))
      if i < #cfl_values then f:write(', ') end
   end
   f:write('],\n')
   f:write('"cfl_schedule_times": [')
   for i,e in ipairs(cfl_times) do
      f:write(string.format('%.18e', e))
      if i < #cfl_times then f:write(', ') end
   end
   f:write('],\n')
   --
   f:write(string.format('"solid_domain_cfl" : %.18e,\n', config.solid_domain_cfl))
   f:write(string.format('"coupling_with_solid_domains": "%s",\n', config.coupling_with_solid_domains))
   f:write(string.format('"solid_domain_augmented_deriv_avg": %s,\n', tostring(config.solid_domain_augmented_deriv_avg)))
   f:write(string.format('"fluid_solid_bc_use_heat_transfer_coeff": %s,\n', tostring(config.fluid_solid_bc_use_heat_transfer_coeff)))
   f:write(string.format('"solid_has_isotropic_properties": %s,\n', tostring(config.solid_has_isotropic_properties)))
   f:write('"solid_thermal_models": {\n')
   for k,v in pairs(_solidModels) do
      f:write(string.format(' "%s": \n', k))
      f:write(string.format('%s,', v:tojson()))
   end
   f:write('\n "dummy_entry_without_trailing_comma": 0\n') -- no comma on last entry
   f:write('},\n')
   --[[ RJG, 2024-02-13 disabled temporarily
   f:write('"solid_domain_loose_update_options" : {\n')
   f:write(string.format('   "max_newton_iterations" : %d,\n', SolidDomainLooseUpdate.max_newton_iterations))
   f:write(string.format('   "newton_solve_tolerance" : %.18e,\n', SolidDomainLooseUpdate.newton_solve_tolerance))
   f:write(string.format('   "max_gmres_iterations" : %d,\n', SolidDomainLooseUpdate.max_gmres_iterations))
   f:write(string.format('   "max_gmres_restarts" : %d,\n', SolidDomainLooseUpdate.max_gmres_restarts))
   f:write(string.format('   "gmres_solve_tolerance" : %.18e,\n', SolidDomainLooseUpdate.gmres_solve_tolerance))
   f:write(string.format('   "perturbation_size" : %.18e,\n', SolidDomainLooseUpdate.perturbation_size))
   f:write(string.format('   "cfl" : %.18e,\n', SolidDomainLooseUpdate.cfl))
   f:write(string.format('   "solid_time_integration_scheme" : "%s",\n', SolidDomainLooseUpdate.solid_time_integration_scheme))
   f:write(string.format('   "solid_domain_only" : %s,\n', tostring(SolidDomainLooseUpdate.solid_domain_only)))
   f:write(string.format('   "super_time_steps" : %d,\n', SolidDomainLooseUpdate.super_time_steps))
   f:write(string.format('   "implicit_time_integration_mode" : %d,\n', SolidDomainLooseUpdate.implicit_time_integration_mode))
   f:write(string.format('   "frozen_preconditioner_count" : %d,\n', SolidDomainLooseUpdate.frozen_preconditioner_count))
   f:write(string.format('   "preconditioner_fill_in" : %d,\n', SolidDomainLooseUpdate.preconditioner_fill_in))
   f:write(string.format('   "preconditioner_approximation" : %d,\n', SolidDomainLooseUpdate.preconditioner_approximation))
   f:write(string.format('   "use_preconditioner" : %s \n', tostring(SolidDomainLooseUpdate.use_preconditioner)))
   -- Note, also, no comma on last entry in JSON object. (^^^: Look up one line and check!)
   f:write('},\n')
   --]]
   f:write(string.format('"MHD": %s,\n', tostring(config.MHD)))
   f:write(string.format('"MHD_static_field": %s,\n', tostring(config.MHD_static_field)))
   f:write(string.format('"MHD_resistive": %s,\n', tostring(config.MHD_resistive)))
   f:write(string.format('"divergence_cleaning": %s,\n', tostring(config.divergence_cleaning)))
   f:write(string.format('"divB_damping_length": %.18e,\n', config.divB_damping_length))
   f:write(string.format('"electric_field_count": %d,\n', config.electric_field_count))
   f:write(string.format('"solve_electric_field": %s,\n', tostring(config.solve_electric_field)))
   f:write(string.format('"field_conductivity_model": "%s",\n', tostring(config.field_conductivity_model)))
   f:write(string.format('"apply_bcs_in_parallel": %s,\n',
			 tostring(config.apply_bcs_in_parallel)))
   f:write(string.format('"flowstate_limits_max_velocity": %.18e,\n', config.flowstate_limits_max_velocity))
   f:write(string.format('"flowstate_limits_max_tke": %.18e,\n', config.flowstate_limits_max_tke))
   f:write(string.format('"flowstate_limits_min_tke": %.18e,\n', config.flowstate_limits_min_tke))
   f:write(string.format('"flowstate_limits_max_temp": %.18e,\n', config.flowstate_limits_max_temp))
   f:write(string.format('"flowstate_limits_min_temp": %.18e,\n', config.flowstate_limits_min_temp))
   f:write(string.format('"flowstate_limits_min_pressure": %.18e,\n', config.flowstate_limits_min_pressure))
   f:write(string.format('"user_specified_velocities_are_in_non_rotating_frame": %s,\n',
                         tostring(config.user_specified_velocities_are_in_non_rotating_frame)))
   f:write(string.format('"max_invalid_cells": %d,\n', config.max_invalid_cells))
   f:write(string.format('"adjust_invalid_cell_data": %s,\n', tostring(config.adjust_invalid_cell_data)))
   f:write(string.format('"report_invalid_cells": %s,\n', tostring(config.report_invalid_cells)))

   if config.n_ghost_cell_layers < 3 then
      if config.high_order_flux_calculator or (config.interpolation_order == 3) then
         print("Increasing config.n_ghost_cell_layers to 3 (structured-grid).")
         config.n_ghost_cell_layers = 3
      end
   end
   f:write(string.format('"n_ghost_cell_layers": %d,\n', config.n_ghost_cell_layers))
   f:write(string.format('"high_order_flux_calculator": %s,\n', tostring(config.high_order_flux_calculator)))
   f:write(string.format('"flux_calculator": "%s",\n', config.flux_calculator))
   f:write(string.format('"interpolation_order": %d,\n', config.interpolation_order))
   f:write(string.format('"interpolation_delay": %.18e,\n', config.interpolation_delay))
   f:write(string.format('"allow_interpolation_for_sensitivity_matrix": %s,\n',
                         tostring(config.allow_interpolation_for_sensitivity_matrix)))
   f:write(string.format('"suppress_radial_reconstruction_at_xaxis": %s,\n',
                         tostring(config.suppress_radial_reconstruction_at_xaxis)))
   f:write(string.format('"suppress_reconstruction_at_shocks": %s,\n',
                         tostring(config.suppress_reconstruction_at_shocks)))
   f:write(string.format('"suppress_reconstruction_at_boundaries": %s,\n',
			 tostring(config.suppress_reconstruction_at_boundaries)))
   f:write(string.format('"thermo_interpolator": "%s",\n',
			 string.lower(config.thermo_interpolator)))
   f:write(string.format('"allow_reconstruction_for_energy_modes": %s,\n',
			 tostring(config.allow_reconstruction_for_energy_modes)))
   f:write(string.format('"allow_reconstruction_for_species": %s,\n',
			 tostring(config.allow_reconstruction_for_species)))
   f:write(string.format('"allow_reconstruction_for_turbulent_variables": %s,\n',
			 tostring(config.allow_reconstruction_for_turbulent_variables)))
   f:write(string.format('"interpolate_in_local_frame": %s,\n',
			 tostring(config.interpolate_in_local_frame)))
   f:write(string.format('"apply_limiter": %s,\n', tostring(config.apply_limiter)))
   f:write(string.format('"epsilon_van_albada": %s,\n', tostring(config.epsilon_van_albada)))
   f:write(string.format('"extrema_clipping": %s,\n', tostring(config.extrema_clipping)))
   f:write(string.format('"apply_heuristic_pressure_based_limiting": %s,\n', tostring(config.apply_heuristic_pressure_based_limiting)))
   f:write(string.format('"apply_entropy_fix": %s,\n', tostring(config.apply_entropy_fix)))
   f:write(string.format('"enforce_species_density_positivity": %s,\n', tostring(config.enforce_species_density_positivity)))
   f:write(string.format('"scale_species_after_reconstruction": %s,\n', tostring(config.scale_species_after_reconstruction)))
   f:write(string.format('"unstructured_limiter": "%s",\n', config.unstructured_limiter))
   f:write(string.format('"freeze_limiter_on_step": %d,\n', config.freeze_limiter_on_step))
   f:write(string.format('"use_extended_stencil": %s,\n', tostring(config.use_extended_stencil)))
   f:write(string.format('"smooth_limiter_coeff": %.18e,\n', config.smooth_limiter_coeff))
   f:write(string.format('"nsteps_of_chemistry_ramp": %d,\n', config.nsteps_of_chemistry_ramp))
   f:write(string.format('"compression_tolerance": %.18e,\n', config.compression_tolerance))
   f:write(string.format('"shock_detector": "%s",\n', config.shock_detector))
   f:write(string.format('"do_shock_detect": %s,\n', tostring(config.do_shock_detect)))
   f:write(string.format('"damped_outflow": %s,\n', tostring(config.damped_outflow)))
   f:write(string.format('"strict_shock_detector": %s,\n', tostring(config.strict_shock_detector)))
   f:write(string.format('"shear_tolerance": %.18e,\n', config.shear_tolerance))
   f:write(string.format('"shock_detector_smoothing": %d,\n', config.shock_detector_smoothing))
   f:write(string.format('"frozen_shock_detector": %s,\n', tostring(config.frozen_shock_detector)))
   f:write(string.format('"shock_detector_freeze_step": %d,\n', config.shock_detector_freeze_step))
   f:write(string.format('"M_inf": %.18e,\n', config.M_inf))
   f:write(string.format('"artificial_compressibility": %s,\n', tostring(config.artificial_compressibility)))
   f:write(string.format('"ac_alpha": %.18e,\n', config.ac_alpha))
   --
   f:write(string.format('"radiation": %s,\n', tostring(config.radiation)))
   --
   f:write(string.format('"grid_motion": "%s",\n', tostring(config.grid_motion)))
   f:write(string.format('"write_vertex_velocities": %s,\n', tostring(config.write_vertex_velocities)))
   f:write(string.format('"udf_grid_motion_file": "%s",\n', tostring(config.udf_grid_motion_file)))
   --
   f:write(string.format('"shock_fitting_delay": %.18e,\n', config.shock_fitting_delay))
   f:write(string.format('"shock_fitting_allow_flow_reconstruction": %s,\n',
                         tostring(config.shock_fitting_allow_flow_reconstruction)))
   f:write(string.format('"shock_fitting_scale_factor": %.18e,\n', config.shock_fitting_scale_factor))
   f:write(string.format('"shock_fitting_filter_velocity_scale": %.18e,\n',
                         config.shock_fitting_filter_velocity_scale))
   f:write(string.format('"shock_fitting_assume_symmetry_at_first_point": %s,\n',
                         tostring(config.shock_fitting_assume_symmetry_at_first_point)))
   --
   f:write(string.format('"viscous": %s,\n', tostring(config.viscous)))
   f:write(string.format('"use_viscosity_from_cells": %s,\n', tostring(config.use_viscosity_from_cells)))
   f:write(string.format('"spatial_deriv_calc": "%s",\n', config.spatial_deriv_calc))
   f:write(string.format('"spatial_deriv_locn": "%s",\n', config.spatial_deriv_locn))
   f:write(string.format('"include_ghost_cells_in_spatial_deriv_clouds": %s,\n',
			 tostring(config.include_ghost_cells_in_spatial_deriv_clouds)))
   f:write(string.format('"upwind_vertex_gradients": %s,\n',tostring(config.upwind_vertex_gradients)))
   f:write(string.format('"save_convective_gradients": %s,\n',tostring(config.save_convective_gradients)))
   f:write(string.format('"save_viscous_gradients": %s,\n',tostring(config.save_viscous_gradients)))
   f:write(string.format('"save_limiter_values": %s,\n',tostring(config.save_limiter_values)))
   f:write(string.format('"save_residual_values": %s,\n',tostring(config.save_residual_values)))
   f:write(string.format('"save_timestep_values": %s,\n',tostring(config.save_timestep_values)))
   f:write(string.format('"nic_write": %d,\n', config.nic_write))
   f:write(string.format('"njc_write": %d,\n', config.njc_write))
   f:write(string.format('"nkc_write": %d,\n', config.nkc_write))
   f:write(string.format('"viscous_factor": %.18e,\n', config.viscous_factor))
   f:write(string.format('"viscous_factor_increment": %.18e,\n', config.viscous_factor_increment))
   f:write(string.format('"viscous_delay": %.18e,\n', config.viscous_delay))
   f:write(string.format('"shear_stress_relative_limit": %.18e,\n', config.shear_stress_relative_limit))
   f:write(string.format('"apply_shear_stress_relative_limit": %s,\n', tostring(config.apply_shear_stress_relative_limit)))
   f:write(string.format('"mass_diffusion_model": "%s",\n', string.lower(config.mass_diffusion_model)))
   f:write(string.format('"diffusion_coefficient_type": "%s",\n', string.lower(config.diffusion_coefficient_type)))
   f:write(string.format('"lewis_number": %.18e,\n', config.lewis_number))
   --
   f:write(string.format('"turbulence_model": "%s",\n', string.lower(config.turbulence_model)))
   f:write(string.format('"turbulence_prandtl_number": %.18e,\n', config.turbulence_prandtl_number))
   f:write(string.format('"turbulence_schmidt_number": %.18e,\n', config.turbulence_schmidt_number))
   f:write(string.format('"max_mu_t_factor": %.18e,\n', config.max_mu_t_factor))
   f:write(string.format('"transient_mu_t_factor": %.18e,\n', config.transient_mu_t_factor))
   f:write(string.format('"freestream_turbulent_intensity": %.18e,\n', config.freestream_turbulent_intensity))
   --
   f:write(string.format('"udf_source_terms_file": "%s",\n', config.udf_source_terms_file))
   f:write(string.format('"udf_source_terms": %s,\n', tostring(config.udf_source_terms)))
   f:write(string.format('"eval_udf_source_terms_at_each_stage": %s,\n', tostring(config.eval_udf_source_terms_at_each_stage)))
   --
   f:write(string.format('"chemistry_update": "%s",\n', config.chemistry_update))
   f:write(string.format('"reacting": %s,\n', tostring(config.reacting)))
   f:write(string.format('"reactions_file": "%s",\n', config.reactions_file))
   f:write(string.format('"reaction_time_delay": %.18e,\n', config.reaction_time_delay))
   if #config.reaction_fraction_schedule > 0 then
      -- The table already have some values.
      -- We will presume that they are valid entries, however,
      -- if enough people get it wrong, we'll put some check here.
   else
      -- Fall back to making up a schedule by assuming full time-step fraction.
      config.reaction_fraction_schedule = {{0.0, 1.0},}
   end
   f:write(string.format('"reaction_fraction_schedule_length": %d,\n', #config.reaction_fraction_schedule))
   local rf_values = {}
   local rf_times = {}
   for i,rf_pair in ipairs(config.reaction_fraction_schedule) do
      rf_times[#rf_times+1] = rf_pair[1]
      rf_values[#rf_values+1] = rf_pair[2]
   end
   f:write('"reaction_fraction_schedule_values": [')
   for i,e in ipairs(rf_values) do
      f:write(string.format('%.18e', e))
      if i < #rf_values then f:write(', ') end
   end
   f:write('],\n')
   f:write('"reaction_fraction_schedule_times": [')
   for i,e in ipairs(rf_times) do
      f:write(string.format('%.18e', e))
      if i < #rf_times then f:write(', ') end
   end
   f:write('],\n')
   f:write(string.format('"T_frozen": %.18e,\n', config.T_frozen))
   f:write(string.format('"T_frozen_energy": %.18e,\n', config.T_frozen_energy))
   f:write(string.format('"tci_model": "%s",\n', string.lower(config.tci_model)))
   f:write(string.format('"ignition_time_start": %.18e,\n', config.ignition_time_start))
   f:write(string.format('"ignition_time_stop": %.18e,\n', config.ignition_time_stop))
   f:write(string.format('"energy_exchange_file": "%s",\n', config.energy_exchange_file))
   f:write(string.format('"radiation_energy_dump_allowed": %s,\n', tostring(config.radiation_energy_dump_allowed)))
   f:write(string.format('"radiation_energy_dump_temperature_limit": %.18e,\n', config.radiation_energy_dump_temperature_limit))
   --
   f:write(string.format('"control_count": %d,\n', config.control_count))
   f:write(string.format('"write_transient_residuals": %s,\n', tostring(config.write_transient_residuals)))
   f:write(string.format('"nfluidblock": %d,\n', #fluidBlocks))
   f:write(string.format('"nfluidblockarrays": %d,\n', #fluidBlockArrays))
   --
   f:write(string.format('"diffuse_wall_bcs_on_init": %s,\n', tostring(config.diffuse_wall_bcs_on_init)))
   f:write(string.format('"number_init_passes": %d,\n', config.number_init_passes))
   f:write(string.format('"wall_temperature_on_init": %.18e,\n', config.wall_temperature_on_init));
   --
   f:write(string.format('"do_temporal_DFT": %s,\n', tostring(config.do_temporal_DFT)))
   f:write(string.format('"DFT_n_modes": %d,\n', config.DFT_n_modes))
   f:write(string.format('"DFT_step_interval": %d,\n', config.DFT_step_interval))
   --
   f:write(string.format('"do_flow_average": %s,\n', tostring(config.do_flow_average)))
   --
   f:write(string.format('"FEMModel": \"%s\",\n', tostring(config.FEMModel)))
   --
   f:write(string.format('"block_marching": %s,\n',
			 tostring(config.block_marching)))
   f:write(string.format('"nib": %d,\n', config.nib))
   f:write(string.format('"njb": %d,\n', config.njb))
   f:write(string.format('"nkb": %d,\n', config.nkb))
   f:write(string.format('"propagate_inflow_data": %s,\n',
			 tostring(config.propagate_inflow_data)))
   f:write(string.format('"save_intermediate_results": %s,\n',
			 tostring(config.save_intermediate_results)))
   f:write(string.format('"boundary_groups_for_loads": "%s",\n',
			 config.boundary_groups_for_loads))
   f:write(string.format('"write_loads": %s,\n',
			 tostring(config.write_loads)))
   f:write(string.format('"compute_run_time_loads": %s,\n',
			 tostring(config.compute_run_time_loads)))
   f:write(string.format('"run_time_loads_count": %d,\n', config.run_time_loads_count))
   if run_time_loads and (#run_time_loads > 0) then
      f:write(string.format('"run_time_loads": { "ngroups" : %d, \n', #run_time_loads))
      for i,t in ipairs(run_time_loads) do
         f:write(string.format('   "group-%d" : { "groupLabel" : "%s",\n',  i-1, run_time_loads[i].group))
         f:write(string.format('       "momentCtr_x" : %.18e,\n',  run_time_loads[i].moment_centre.x))
         f:write(string.format('       "momentCtr_y" : %.18e,\n',  run_time_loads[i].moment_centre.y))
         f:write(string.format('       "momentCtr_z" : %.18e },\n',  run_time_loads[i].moment_centre.z))
      end
      f:write('   "dummy_entry_without_trailing_comma": 0\n') -- no comma on last entry
      f:write('},\n')
   end
   f:write(string.format('"nhcell": %d,\n', #historyCells))
   for i,hcell in ipairs(historyCells) do
      f:write(string.format('"history-cell-%d": [%d, %d],\n', i-1, hcell.ib, hcell.i))
   end
   f:write(string.format('"nsolidhcell": %d,\n', #solidHistoryCells))
   for i,hcell in ipairs(solidHistoryCells) do
      f:write(string.format('"solid-history-cell-%d": [%d, %d],\n', i-1, hcell.ib, hcell.i))
   end
   --
   f:write(string.format('"n-reaction-zones": %d,\n', #reactionZones))
   for i,zone in ipairs(reactionZones) do
      f:write(string.format('"reaction-zone-%d": [%.18e, %.18e, %.18e, %.18e, %.18e, %.18e],\n',
			    i-1, zone.p0.x, zone.p0.y, zone.p0.z,
			    zone.p1.x, zone.p1.y, zone.p1.z))
   end
   f:write(string.format('"n-ignition-zones": %d,\n', #ignitionZones))
   for i,zone in ipairs(ignitionZones) do
      f:write(string.format('"ignition-zone-%d": [%.18e, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e],\n',
			    i-1, zone.p0.x, zone.p0.y, zone.p0.z,
			    zone.p1.x, zone.p1.y, zone.p1.z, zone.T))
   end
   f:write(string.format('"n-turbulent-zones": %d,\n', #turbulentZones))
   for i,zone in ipairs(turbulentZones) do
      f:write(string.format('"turbulent-zone-%d": [%.18e, %.18e, %.18e, %.18e, %.18e, %.18e],\n',
			    i-1, zone.p0.x, zone.p0.y, zone.p0.z,
			    zone.p1.x, zone.p1.y, zone.p1.z))
   end
   f:write(string.format('"n-suppress-reconstruction-zones": %d,\n', #suppressReconstructionZones))
   for i,zone in ipairs(suppressReconstructionZones) do
      f:write(string.format('"suppress-reconstruction-zone-%d": [%.18e, %.18e, %.18e, %.18e, %.18e, %.18e],\n',
			    i-1, zone.p0.x, zone.p0.y, zone.p0.z,
			    zone.p1.x, zone.p1.y, zone.p1.z))
   end
   f:write(string.format('"n-suppress-viscous-stresses-zones": %d,\n', #suppressViscousStressesZones))
   for i,zone in ipairs(suppressViscousStressesZones) do
      f:write(string.format('"suppress-viscous-stresses-zone-%d": [%.18e, %.18e, %.18e, %.18e, %.18e, %.18e],\n',
			    i-1, zone.p0.x, zone.p0.y, zone.p0.z,
			    zone.p1.x, zone.p1.y, zone.p1.z))
   end
   --
   f:write(string.format('"udf_solid_source_terms_file": "%s",\n', config.udf_solid_source_terms_file))
   f:write(string.format('"udf_solid_source_terms": %s,\n', tostring(config.udf_solid_source_terms)))
   f:write(string.format('"nsolidblock": %d,\n', #solidBlocks))
   --
   for i = 1, #fluidBlockArrays do
      f:write(string.format('"fluid_block_array_%d": ', (i-1)))
      f:write(json.stringify(fluidBlockArrays[i]) .. ",\n")
   end
   for i = 1, #fluidBlocks do
      f:write(fluidBlocks[i]:tojson() .. ",\n")
   end
   for i = 1, #solidBlocks do
      f:write(solidBlocks[i]:tojson() .. ",\n")
   end
   --
   f:write('"dummy_entry_without_trailing_comma": 0\n') -- no comma on last entry
   f:write("}\n")
   --
   f:close()
end

function output.write_times_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# tindx sim_time dt_global\n");
   f:write(string.format("%04d %.18e %.18e\n", 0, config.start_time, config.dt_init))
   f:close()
end

function output.write_block_list_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# indx type label ncells\n")
   for i = 1, #(fluidBlocks) do
      local blk = fluidBlocks[i]
      local grid_type
      if blk.grid then
         grid_type = blk.grid:get_type()
      else
         grid_type = blk.gridMetadata.type
      end
      f:write(string.format("%4d %s %s %d\n", blk.id, grid_type, blk.label, blk.ncells))
   end
   f:close()
end

function output.write_mpimap_file(fileName)
   if not mpiTasks then
      -- The user's input script has not set up mpiTasks, so we need to do it now.
      if config.block_marching then
         -- Work through fluidBlockArrays and allocate MPI tasks for each block array.
         for _,fba in ipairs(fluidBlockArrays) do
            mpiDistributeFluidBlockArray{fba=fba, ntasks=fba.njb*fba.nkb}
         end
      else
         -- Work through the single-dimensional fluidBlocks list
         mpiDistributeBlocks()
      end
   end
   local f = assert(io.open(fileName, "w"))
   f:write("# indx mpiTask\n")
   for i = 1, #(fluidBlocks) do
      f:write(string.format("%4d %4d\n", fluidBlocks[i].id, mpiTasks[i]))
   end
   for i = 1, #(solidBlocks) do
      f:write(string.format("%4d %4d\n", solidBlocks[i].id, mpiTasks[i + #fluidBlocks]))
   end
   f:close()
end

return output
