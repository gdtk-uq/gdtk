
LMR ?= .
LMR_CMD = $(LMR)/commands
LMR_LUA_MOD = $(LMR)/lua-modules
LMR_LUA_WRAP = $(LMR)/luawrap

LMR_CORE_FILES = $(LMR)/block.d \
	$(LMR)/coredata.d \
	$(LMR)/blockio.d \
	$(LMR)/conservedquantities.d \
	$(LMR)/fileutil.d \
	$(LMR)/flowgradients.d \
	$(LMR)/flowsolution.d \
	$(LMR)/flowstate.d \
	$(LMR)/fluidblock.d \
	$(LMR)/fluidblockarray.d \
	$(LMR)/fluxcalc.d \
	$(LMR)/fluidfvcell.d \
	$(LMR)/fvcell.d \
	$(LMR)/fvcellio.d \
	$(LMR)/fvinterface.d \
	$(LMR)/fvvertex.d \
	$(LMR)/globalconfig.d \
	$(LMR)/globaldata.d \
	$(LMR)/grid_motion.d \
	$(LMR)/grid_motion_shock_fitting.d \
	$(LMR)/grid_motion_udf.d \
	$(LMR)/history.d \
	$(LMR)/init.d \
	$(LMR)/jacobian.d \
	$(LMR)/lmrerrors.d \
	$(LMR)/lmrexceptions.d \
	$(LMR)/loads.d \
	$(LMR)/lsqinterp.d \
	$(LMR)/lua_helper.d \
	$(LMR)/mass_diffusion.d \
	$(LMR)/newtonkrylovsolver.d \
	$(LMR)/onedinterp.d \
	$(LMR)/sfluidblock.d \
	$(LMR)/shockdetectors.d \
	$(LMR)/simcore.d \
	$(LMR)/simcore_exchange.d \
	$(LMR)/simcore_gasdynamic_step.d \
	$(LMR)/simcore_solid_step.d \
	$(LMR)/special_block_init.d \
	$(LMR)/timemarching.d \
	$(LMR)/turbulence.d \
	$(LMR)/ufluidblock.d \
	$(LMR)/user_defined_source_terms.d \
	$(LMR)/vtk_writer.d

LMR_BC_FILES = $(LMR)/bc/package.d \
	$(LMR)/bc/boundary_condition.d \
	$(LMR)/bc/ghost_cell_effect/package.d \
	$(LMR)/bc/ghost_cell_effect/ghost_cell.d \
	$(LMR)/bc/ghost_cell_effect/internal_copy_then_reflect.d \
	$(LMR)/bc/ghost_cell_effect/flow_state_copy.d \
	$(LMR)/bc/ghost_cell_effect/flow_state_copy_from_profile.d \
	$(LMR)/bc/ghost_cell_effect/flow_state_copy_from_history.d \
	$(LMR)/bc/ghost_cell_effect/synthesise_flow_state.d \
	$(LMR)/bc/ghost_cell_effect/extrapolate_copy.d \
	$(LMR)/bc/ghost_cell_effect/from_upwind.d \
	$(LMR)/bc/ghost_cell_effect/fixed_p.d \
	$(LMR)/bc/ghost_cell_effect/fixed_pt.d \
	$(LMR)/bc/ghost_cell_effect/from_stagnation.d \
	$(LMR)/bc/ghost_cell_effect/full_face_copy.d \
	$(LMR)/bc/ghost_cell_effect/mapped_cell_copy.d \
	$(LMR)/bc/ghost_cell_effect/gas_solid_full_face_copy.d \
	$(LMR)/bc/user_defined_effects.d \
	$(LMR)/bc/boundary_flux_effect.d \
	$(LMR)/bc/boundary_cell_effect.d \
	$(LMR)/bc/boundary_interface_effect.d

LMR_SOLID_FILES := $(LMR)/solid/solidbc.d \
	$(LMR)/solid/solidblock.d \
	$(LMR)/solid/solid_ghost_cell.d \
	$(LMR)/solid/solid_boundary_flux_effect.d \
	$(LMR)/solid/solid_boundary_interface_effect.d \
	$(LMR)/solid/solid_gas_full_face_copy.d \
	$(LMR)/solid/solid_full_face_copy.d \
	$(LMR)/solid/ssolidblock.d \
	$(LMR)/solid/solidfvcell.d \
	$(LMR)/solid/solidfvinterface.d \
	$(LMR)/solid/solidfvvertex.d \
	$(LMR)/solid/solidsolution.d \
	$(LMR)/solid/solidstate.d \
	$(LMR)/solid/solidthermalmodel.d \
	$(LMR)/solid/solid_udf_source_terms.d

LMR_EFIELD_FILES := $(LMR)/efield/efield.d \
	$(LMR)/efield/efieldgmres.d \
	$(LMR)/efield/efieldconductivity.d \
	$(LMR)/efield/efieldexchange.d \
	$(LMR)/efield/efieldderivatives.d \
	$(LMR)/efield/efieldbc.d

LMR_LUA_FILES = $(LMR_LUA_WRAP)/luaflowsolution.d \
	$(LMR_LUA_WRAP)/luaflowstate.d

LMR_CMD_FILES = $(LMR_CMD)/cmdhelper.d \
	$(LMR_CMD)/command.d \
	$(LMR_CMD)/computenorms.d \
	$(LMR_CMD)/customscript.d \
	$(LMR_CMD)/probeflow.d \
	$(LMR_CMD)/sliceflow.d \
	$(LMR_CMD)/slicesolid.d \
	$(LMR_CMD)/extractline.d \
	$(LMR_CMD)/listspecies.d \
	$(LMR_CMD)/limiter2vtk.d \
	$(LMR_CMD)/residual2vtk.d \
	$(LMR_CMD)/plotdiagnostics.d \
	$(LMR_CMD)/prepenergyexchange.d \
	$(LMR_CMD)/prepgas.d \
	$(LMR_CMD)/prepsim.d \
	$(LMR_CMD)/prepgrids.d \
	$(LMR_CMD)/prepmappedcells.d \
	$(LMR_CMD)/prepreactions.d \
	$(LMR_CMD)/revisionid.d \
	$(LMR_CMD)/runsim.d \
	$(LMR_CMD)/snapshot2vtk.d \
	$(LMR_CMD)/structured2unstructured.d
