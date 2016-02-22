KINETICS_DIR ?= .
KINETICS_FILES := $(KINETICS_DIR)/package.d \
	$(KINETICS_DIR)/chemistry_update.d \
	$(KINETICS_DIR)/rate_constant.d \
	$(KINETICS_DIR)/reaction.d \
	$(KINETICS_DIR)/reaction_mechanism.d

KINETICS_LUA_FILES := $(KINETICS_DIR)/luachemistry_update.d \
	$(KINETICS_DIR)/luareaction_mechanism.d
