
LMR ?= .
LMR_CMD = $(LMR)/commands
LMR_LUA_MOD = $(LMR)/lua-modules

LMR_CORE_FILES = $(LMR)/flowsolution.d \
	$(LMR)/fluidblock.d \
	$(LMR)/fvcell.d \
	$(LMR)/jacobian.d \
	$(LMR)/lmrconfig.d \
	$(LMR)/newtonkrylovsolver.d

LMR_CMD_FILES = $(LMR_CMD)/command.d \
	$(LMR_CMD)/prepflow.d \
	$(LMR_CMD)/prepgrids.d \
	$(LMR_CMD)/runsteady.d \
	$(LMR_CMD)/snapshot2vtk.d

LMR_LUA_MODULES = $(LMR_LUA_MOD)/lmrconfig.lua \
	$(LMR_LUA_MOD)/nkconfig.lua
