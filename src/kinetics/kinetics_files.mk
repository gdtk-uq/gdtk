KINETICS_DIR ?= .
KINETICS_FILES := $(KINETICS_DIR)/package.d \
	$(KINETICS_DIR)/thermochemical_reactor.d \
	$(KINETICS_DIR)/init_thermochemical_reactor.d \
	$(KINETICS_DIR)/chemistry_update.d \
	$(KINETICS_DIR)/electronic_update.d \
	$(KINETICS_DIR)/electronically_specific_kinetics.d \
	$(KINETICS_DIR)/ideal_dissociating_gas_kinetics.d \
	$(KINETICS_DIR)/fuel_air_mix_kinetics.d \
	$(KINETICS_DIR)/powers_aslam_kinetics.d \
	$(KINETICS_DIR)/rate_constant.d \
	$(KINETICS_DIR)/reaction.d \
	$(KINETICS_DIR)/reaction_mechanism.d \
	$(KINETICS_DIR)/two_temperature_air_kinetics.d \
	$(KINETICS_DIR)/two_temperature_argon_kinetics.d \
	$(KINETICS_DIR)/two_temperature_nitrogen_kinetics.d \
	$(KINETICS_DIR)/vib_specific_nitrogen_kinetics.d \
	$(KINETICS_DIR)/two_temperature_gasgiant_kinetics.d

ifeq ($(WITH_DVODE),1)
    KINETICS_FILES += $(KINETICS_DIR)/pseudo_species_kinetics.d \
	$(KINETICS_DIR)/dvode.o \
	$(KINETICS_DIR)/pseudosp_solve_ode.o \
	$(KINETICS_DIR)/pseudosp_rhs_ode.o 
endif

KINETICS_LUA_FILES := $(KINETICS_DIR)/luathermochemical_reactor.d \
	$(KINETICS_DIR)/luachemistry_update.d \
	$(KINETICS_DIR)/luaelectronically_specific_kinetics.d \
	$(KINETICS_DIR)/luareaction_mechanism.d \
	$(KINETICS_DIR)/luatwo_temperature_air_kinetics.d \
	$(KINETICS_DIR)/luavib_specific_nitrogen_kinetics.d

ifeq ($(WITH_DVODE),1)
    KINETICS_LUA_FILES += $(KINETICS_DIR)/luapseudo_species_kinetics.d
endif
