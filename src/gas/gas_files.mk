GAS_DIR ?= .
GAS_MODEL_FILES := $(GAS_DIR)/package.d \
	$(GAS_DIR)/co2gas.d \
	$(GAS_DIR)/co2gas_sw.d \
	$(GAS_DIR)/gas_model.d \
	$(GAS_DIR)/gas_state.d \
	$(GAS_DIR)/ideal_gas.d \
	$(GAS_DIR)/cea_gas.d \
	$(GAS_DIR)/physical_constants.d \
	$(GAS_DIR)/sf6virial.d \
	$(GAS_DIR)/therm_perf_gas.d \
	$(GAS_DIR)/very_viscous_air.d \
	$(GAS_DIR)/uniform_lut.d \
	$(GAS_DIR)/adaptive_lut_CEA.d \
	$(GAS_DIR)/ideal_air_proxy.d \
	$(GAS_DIR)/powers_aslam_gas.d \
	$(GAS_DIR)/two_temperature_reacting_argon.d \
	$(GAS_DIR)/ideal_dissociating_gas.d \
	$(GAS_DIR)/two_temperature_nitrogen.d \
	$(GAS_DIR)/vib_specific_nitrogen.d \
	$(GAS_DIR)/fuel_air_mix.d \
	$(GAS_DIR)/ideal_air_fortran.o \
	$(GAS_DIR)/steam.d

THERMO_FILES := \
	$(GAS_DIR)/thermo/cea_thermo_curves.d \
	$(GAS_DIR)/thermo/evt_eos.d \
	$(GAS_DIR)/thermo/perf_gas_mix_eos.d \
	$(GAS_DIR)/thermo/pvt_eos.d \
	$(GAS_DIR)/thermo/therm_perf_gas_mix_eos.d

DIFFUSION_FILES := \
	$(GAS_DIR)/diffusion/cea_therm_cond.d \
	$(GAS_DIR)/diffusion/cea_viscosity.d \
	$(GAS_DIR)/diffusion/sutherland_therm_cond.d \
	$(GAS_DIR)/diffusion/sutherland_viscosity.d \
	$(GAS_DIR)/diffusion/therm_cond.d \
	$(GAS_DIR)/diffusion/viscosity.d \
	$(GAS_DIR)/diffusion/wilke_mixing_therm_cond.d \
	$(GAS_DIR)/diffusion/wilke_mixing_viscosity.d

GAS_FILES := $(GAS_MODEL_FILES) $(THERMO_FILES) $(DIFFUSION_FILES)

GAS_LUA_FILES := $(GAS_DIR)/luagas_model.d
