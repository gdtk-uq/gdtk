# makefile for the diffusion module
UTIL_DIR := ../../util
UTIL_FILES := $(UTIL_DIR)/lua.d \
	$(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

LUA := ../../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib
LUAD_DIR := ../../extern/LuaD
LUAD_FILES := $(LUAD_DIR)/luad/*.d \
	$(LUAD_DIR)/luad/c/*.d \
	$(LUAD_DIR)/luad/conversions/*.d

GAS_DIR := ..
LIBGAS := ${GAS_DIR}/libgas.a

DLINKFLAGS :=  -w -L-L${LIBLUAPATH} -L-llua -L-ldl

unittest: sutherland_visc sutherland_therm_cond \
	wilke_mixing_visc wilke_mixing_therm_cond \
	cea_viscosity cea_therm_cond
	./sutherland_visc
	./sutherland_therm_cond
	./wilke_mixing_visc
	./wilke_mixing_therm_cond
	./cea_viscosity

${LIBGAS}:
	cd ${GAS_DIR}; make libgas.a

${LIBLUA}:
	cd ${LUA_DIR}; make linux local

sutherland_viscosity : sutherland_viscosity.d viscosity.d ../gas_model.d $(UTIL_FILES) $(LIBLUA)
	dmd -main -unittest sutherland_viscosity.d viscosity.d ../gas_model.d ../physical_constants.d -I../ $(UTIL_FILES) -L-L$(LIBLUAPATH) -L-llua -L-ldl

sutherland_therm_cond : sutherland_therm_cond.d therm_cond.d ../gas_model.d $(UTIL_FILES) $(LIBLUA)
	dmd -main -unittest sutherland_therm_cond.d therm_cond.d ../gas_model.d ../physical_constants.d -I../ $(UTIL_FILES) -L-L$(LIBLUAPATH) -L-llua -L-ldl

wilke_mixing_viscosity : wilke_mixing_viscosity.d $(LIBGAS) $(LIBLUA)
	dmd -main -unittest wilke_mixing_viscosity.d $(LIBGAS) $(LUAD_FILES) $(UTIL_FILES) -I../../ -I$(UTIL_DIR) $(DLINKFLAGS)

wilke_mixing_therm_cond : wilke_mixing_therm_cond.d $(LIBGAS) $(LIBLUA)
	dmd -main -unittest wilke_mixing_therm_cond.d $(LIBGAS) $(LUAD_FILES) $(UTIL_FILES) -I../../ -I$(UTIL_DIR) $(DLINKFLAGS)

cea_viscosity : cea_viscosity.d viscosity.d ../gasmodel.d $(UTIL_FILES) $(LIBLUA)
	dmd -main -unittest cea_viscosity.d viscosity.d ../gasmodel.d -I../ $(UTIL_FILES) ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

cea_therm_cond : cea_therm_cond.d therm_cond.d ../gasmodel.d $(UTIL_FILES) $(LIBLUA)
	dmd -main -unittest cea_therm_cond.d therm_cond.d ../gasmodel.d -I../ $(UTIL_FILES) ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

