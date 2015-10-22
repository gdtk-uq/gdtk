GAS_MODEL_FILES := package.d \
	co2gas.d \
	co2gas_sw.d \
	gas_model.d \
	gas_model_util.d \
	ideal_gas.d \
	physical_constants.d \
	sf6virial.d \
	therm_perf_gas.d \
	very_viscous_air.d


THERMO_FILES := \
	thermo/cea_thermo_curves.d \
	thermo/evt_eos.d \
	thermo/perf_gas_mix_eos.d \
	thermo/pvt_eos.d \
	thermo/therm_perf_gas_mix_eos.d

DIFFUSION_FILES := \
	diffusion/cea_therm_cond.d \
	diffusion/cea_viscosity.d \
	diffusion/sutherland_therm_cond.d \
	diffusion/sutherland_viscosity.d \
	diffusion/therm_cond.d \
	diffusion/viscosity.d \
	diffusion/wilke_mixing_therm_cond.d \
	diffusion/wilke_mixing_viscosity.d

GAS_FILES := $(GAS_MODEL_FILES) $(THERMO_FILES) $(DIFFUSION_FILES)
