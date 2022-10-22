LMR4 ?= ../eilmer

# Try to add to list in alphabetical order.
LMR4_CORE_FILES := $(LMR4)/block.d \
	$(LMR4)/celldata.d \
	$(LMR4)/conservedquantities.d \
	$(LMR4)/fileutil.d \
	$(LMR4)/flowgradients.d \
	$(LMR4)/flowstate.d \
	$(LMR4)/fluidblockarray.d \
	$(LMR4)/fluidblockio.d \
	$(LMR4)/fluidblockio_new.d \
	$(LMR4)/fluidblockio_old.d \
	$(LMR4)/fluxcalc.d \
	$(LMR4)/fvinterface.d \
	$(LMR4)/fvvertex.d \
	$(LMR4)/gas_solid_interface.d \
	$(LMR4)/globalconfig.d \
	$(LMR4)/globaldata.d \
	$(LMR4)/grid_motion.d \
	$(LMR4)/grid_motion_shock_fitting.d \
	$(LMR4)/grid_motion_udf.d \
	$(LMR4)/history.d \
	$(LMR4)/loads.d \
	$(LMR4)/lsqinterp.d \
	$(LMR4)/mass_diffusion.d \
	$(LMR4)/onedinterp.d \
	$(LMR4)/sfluidblock.d \
	$(LMR4)/shockdetectors.d \
	$(LMR4)/simcore.d \
	$(LMR4)/simcore_exchange.d \
	$(LMR4)/simcore_gasdynamic_step.d \
	$(LMR4)/simcore_io.d \
	$(LMR4)/simcore_solid_step.d \
	$(LMR4)/special_block_init.d \
	$(LMR4)/turbulence.d \
	$(LMR4)/ufluidblock.d

LMR4_LUA_FILES := $(LMR4)/luaflowsolution.d \
	$(LMR4)/luaflowstate.d \
	$(LMR4)/lua_helper.d \
	$(LMR4)/user_defined_source_terms.d

LMR4_BC_FILES := $(LMR4)/bc/package.d \
	$(LMR4)/bc/boundary_condition.d \
	$(LMR4)/bc/ghost_cell_effect/package.d \
	$(LMR4)/bc/ghost_cell_effect/ghost_cell.d \
	$(LMR4)/bc/ghost_cell_effect/internal_copy_then_reflect.d \
	$(LMR4)/bc/ghost_cell_effect/flow_state_copy.d \
	$(LMR4)/bc/ghost_cell_effect/flow_state_copy_from_profile.d \
	$(LMR4)/bc/ghost_cell_effect/flow_state_copy_from_history.d \
	$(LMR4)/bc/ghost_cell_effect/synthesise_flow_state.d \
	$(LMR4)/bc/ghost_cell_effect/extrapolate_copy.d \
	$(LMR4)/bc/ghost_cell_effect/from_upwind.d \
	$(LMR4)/bc/ghost_cell_effect/fixed_p.d \
	$(LMR4)/bc/ghost_cell_effect/fixed_pt.d \
	$(LMR4)/bc/ghost_cell_effect/from_stagnation.d \
	$(LMR4)/bc/ghost_cell_effect/full_face_copy.d \
	$(LMR4)/bc/ghost_cell_effect/mapped_cell_copy.d \
	$(LMR4)/bc/ghost_cell_effect/gas_solid_full_face_copy.d \
	$(LMR4)/bc/user_defined_effects.d \
	$(LMR4)/bc/boundary_flux_effect.d \
	$(LMR4)/bc/boundary_cell_effect.d \
	$(LMR4)/bc/boundary_interface_effect.d

LMR4_SOLID_FILES := $(LMR4)/solid/solidbc.d \
	$(LMR4)/solid/solidblock.d \
	$(LMR4)/solid/solid_ghost_cell.d \
	$(LMR4)/solid/solid_boundary_flux_effect.d \
	$(LMR4)/solid/solid_boundary_interface_effect.d \
	$(LMR4)/solid/solid_gas_full_face_copy.d \
	$(LMR4)/solid/solid_full_face_copy.d \
	$(LMR4)/solid/ssolidblock.d \
	$(LMR4)/solid/solidfvcell.d \
	$(LMR4)/solid/solidfvinterface.d \
	$(LMR4)/solid/solidfvvertex.d \
	$(LMR4)/solid/solidprops.d \
	$(LMR4)/solid/solidsolution.d \
	$(LMR4)/solid/solid_udf_source_terms.d \
	$(LMR4)/solid/luasolidprops.d

LMR4_LUA_MODULES := $(LMR4)/configoptions.lua \
	$(LMR4)/blk_conn.lua \
	$(LMR4)/bc.lua \
	$(LMR4)/gridpro.lua \
	$(LMR4)/grid.lua \
	$(LMR4)/gridarray.lua \
	$(LMR4)/flowstate.lua \
	$(LMR4)/fluidblock.lua \
	$(LMR4)/fbarray.lua \
	$(LMR4)/lmr_config.lua \
	$(LMR4)/solidblock.lua \
	$(LMR4)/mpi.lua \
	$(LMR4)/history.lua \
	$(LMR4)/zones.lua \
	$(LMR4)/output.lua \
	$(LMR4)/sssoptions.lua \
	$(LMR4)/prep_check.lua


LMR4_EFIELD_FILES := $(LMR4)/field/field.d \
	$(LMR4)/field/fieldgmres.d \
	$(LMR4)/field/fieldconductivity.d \
	$(LMR4)/field/fieldexchange.d \
	$(LMR4)/field/fieldbc.d

LMR4_EXTRA_FILES := $(LMR4)/postprocess.d \
	$(LMR4)/vtk_writer.d \
	$(LMR4)/tecplot_writer_classic.d \
	$(LMR4)/tecplot_writer.d	
