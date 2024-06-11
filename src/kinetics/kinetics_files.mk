KINETICS_DIR ?= .
KINETICS_FILES := $(KINETICS_DIR)/package.d \
	$(KINETICS_DIR)/thermochemical_reactor.d \
	$(KINETICS_DIR)/init_thermochemical_reactor.d \
	$(KINETICS_DIR)/chemistry_update.d \
	$(KINETICS_DIR)/energy_exchange_mechanism.d \
	$(KINETICS_DIR)/energy_exchange_system.d \
	$(KINETICS_DIR)/equilibrium_update.d \
	$(KINETICS_DIR)/ideal_dissociating_gas_kinetics.d \
	$(KINETICS_DIR)/fuel_air_mix_kinetics.d \
	$(KINETICS_DIR)/powers_aslam_kinetics.d \
	$(KINETICS_DIR)/yee_kotov_kinetics.d \
	$(KINETICS_DIR)/rate_constant.d \
	$(KINETICS_DIR)/reaction.d \
	$(KINETICS_DIR)/reaction_mechanism.d \
	$(KINETICS_DIR)/relaxation_time.d \
	$(KINETICS_DIR)/exchange_cross_section.d \
	$(KINETICS_DIR)/exchange_chemistry_coupling.d \
	$(KINETICS_DIR)/multi_temperature_thermochemical_reactor.d \
	$(KINETICS_DIR)/two_temperature_air_kinetics.d \
	$(KINETICS_DIR)/two_temperature_argon_kinetics.d \
	$(KINETICS_DIR)/two_temperature_argon_with_ideal_gas.d \
	$(KINETICS_DIR)/two_temperature_nitrogen_kinetics.d \
	$(KINETICS_DIR)/two_temperature_dissociating_nitrogen_kinetics.d \
	$(KINETICS_DIR)/vib_specific_nitrogen_kinetics.d \
	$(KINETICS_DIR)/vib_specific_co_kinetics.d \
	$(KINETICS_DIR)/two_temperature_gasgiant_kinetics.d

KINETICS_LUA_FILES := $(KINETICS_DIR)/luathermochemical_reactor.d \
	$(KINETICS_DIR)/luachemistry_update.d \
	$(KINETICS_DIR)/luaequilibrium_calculator.d \
	$(KINETICS_DIR)/luareaction_mechanism.d \
	$(KINETICS_DIR)/luatwo_temperature_air_kinetics.d \
	$(KINETICS_DIR)/luavib_specific_nitrogen_kinetics.d
