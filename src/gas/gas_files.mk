# gas_files.mk
#
# We list the source files that are needed to build the gas models package.
#
GAS_DIR ?= .
GAS_MODEL_FILES := $(GAS_DIR)/package.d \
	$(GAS_DIR)/composite_gas.d \
	$(GAS_DIR)/gas_model.d \
	$(GAS_DIR)/gas_state.d \
	$(GAS_DIR)/init_gas_model.d \
	$(GAS_DIR)/ideal_gas.d \
	$(GAS_DIR)/ideal_helium.d \
	$(GAS_DIR)/cubic_gas.d \
	$(GAS_DIR)/cea_gas.d \
	$(GAS_DIR)/physical_constants.d \
	$(GAS_DIR)/therm_perf_gas.d \
	$(GAS_DIR)/therm_perf_gas_equil.d \
	$(GAS_DIR)/very_viscous_air.d \
	$(GAS_DIR)/uniform_lut.d \
	$(GAS_DIR)/uniform_lut_plus_ideal.d \
	$(GAS_DIR)/adaptive_lut_CEA.d \
	$(GAS_DIR)/ideal_gas_ab.d \
	$(GAS_DIR)/two_temperature_reacting_argon.d \
	$(GAS_DIR)/two_temperature_argon_plus_ideal.d \
	$(GAS_DIR)/ideal_dissociating_gas.d \
	$(GAS_DIR)/two_temperature_air.d \
	$(GAS_DIR)/two_temperature_nitrogen.d \
	$(GAS_DIR)/two_temperature_dissociating_nitrogen.d \
	$(GAS_DIR)/vib_specific_nitrogen.d \
	$(GAS_DIR)/vib_specific_co.d \
	$(GAS_DIR)/fuel_air_mix.d \
	$(GAS_DIR)/equilibrium_gas.d \
	$(GAS_DIR)/two_temperature_gasgiant.d

THERMO_FILES := $(GAS_DIR)/thermo/package.d \
	$(GAS_DIR)/thermo/cea_thermo_curves.d \
	$(GAS_DIR)/thermo/evt_eos.d \
	$(GAS_DIR)/thermo/perf_gas_mix_eos.d \
	$(GAS_DIR)/thermo/pvt_eos.d \
	$(GAS_DIR)/thermo/therm_perf_gas_mix_eos.d \
	$(GAS_DIR)/thermo/thermo_model.d \
	$(GAS_DIR)/thermo/therm_perf_gas_mix.d \
	$(GAS_DIR)/thermo/two_temperature_gas.d \
	$(GAS_DIR)/thermo/three_temperature_gas.d \
	$(GAS_DIR)/thermo/multi_temperature_gas.d \
	$(GAS_DIR)/thermo/energy_modes.d

DIFFUSION_FILES := $(GAS_DIR)/diffusion/package.d \
	$(GAS_DIR)/diffusion/cea_therm_cond.d \
	$(GAS_DIR)/diffusion/cea_viscosity.d \
	$(GAS_DIR)/diffusion/chemkin_therm_cond.d \
	$(GAS_DIR)/diffusion/chemkin_viscosity.d \
	$(GAS_DIR)/diffusion/gas_mixtures.d \
	$(GAS_DIR)/diffusion/sutherland_therm_cond.d \
	$(GAS_DIR)/diffusion/sutherland_viscosity.d \
	$(GAS_DIR)/diffusion/therm_cond.d \
	$(GAS_DIR)/diffusion/transport_properties_model.d \
	$(GAS_DIR)/diffusion/two_temperature_trans_props.d \
	$(GAS_DIR)/diffusion/multi_temperature_trans_props.d \
	$(GAS_DIR)/diffusion/three_temperature_trans_props.d \
	$(GAS_DIR)/diffusion/viscosity.d \
	$(GAS_DIR)/diffusion/wilke_mixing_therm_cond.d \
	$(GAS_DIR)/diffusion/wilke_mixing_viscosity.d \
	$(GAS_DIR)/diffusion/gasgiant_transport_properties.d \
	$(GAS_DIR)/diffusion/binary_diffusion_coefficients.d \
	$(GAS_DIR)/diffusion/rps_diffusion_coefficients.d

GAS_FILES := $(GAS_MODEL_FILES) $(THERMO_FILES) $(DIFFUSION_FILES)

GAS_LUA_FILES := $(GAS_DIR)/luagas_model.d
