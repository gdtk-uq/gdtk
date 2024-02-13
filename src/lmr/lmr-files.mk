
LMR ?= .
LMR_CMD = $(LMR)/commands
LMR_LUA_MOD = $(LMR)/lua-modules
LMR_LUA_WRAP = $(LMR)/luawrap

LMR_CORE_FILES = $(LMR)/blockio.d \
	$(LMR)/flowsolution.d \
	$(LMR)/flowstate.d \
	$(LMR)/fluidblock.d \
	$(LMR)/fvcell.d \
	$(LMR)/fvcellio.d \
	$(LMR)/globalconfig.d \
	$(LMR)/globaldata.d \
	$(LMR)/history.d \
	$(LMR)/init.d \
	$(LMR)/jacobian.d \
	$(LMR)/lmrexceptions.d \
	$(LMR)/loads.d \
	$(LMR)/newtonkrylovsolver.d \
	$(LMR)/sfluidblock.d \
	$(LMR)/simcore.d \
	$(LMR)/timemarching.d \
	$(LMR)/ufluidblock.d

LMR_LUA_FILES = $(LMR_LUA_WRAP)/luaflowsolution.d \
	$(LMR_LUA_WRAP)/luaflowstate.d

LMR_CMD_FILES = $(LMR_CMD)/checkjacobian.d \
	$(LMR_CMD)/cmdhelper.d \
	$(LMR_CMD)/command.d \
	$(LMR_CMD)/computenorms.d \
	$(LMR_CMD)/limiter2vtk.d \
	$(LMR_CMD)/prepsim.d \
	$(LMR_CMD)/prepgrids.d \
	$(LMR_CMD)/prepmappedcells.d \
	$(LMR_CMD)/revisionid.d \
	$(LMR_CMD)/runsim.d \
	$(LMR_CMD)/snapshot2vtk.d \
	$(LMR_CMD)/structured2unstructured.d
