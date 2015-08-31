
UTIL_DIR := ../../util
UTIL_FILES := $(UTIL_DIR)/lua.d \
	$(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

NM_DIR := ../../nm
NM_FILES := $(NM_DIR)/ridder.d

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

unittest: cea_thermo_curves therm_perf_gas_mix_EOS
	./cea_thermo_curves
	./therm_perf_gas_mix_EOS

${LIBGAS}:
	cd ${GAS_DIR}; make libgas.a

${LIBLUA}:
	cd ${LUA_DIR}; make linux local

cea_thermo_curves : cea_thermo_curves.d $(LIBLUA)
	dmd -main -unittest cea_thermo_curves.d $(UTIL_FILES) -I../../ -I$(UTIL_DIR) $(DLINKFLAGS)

therm_perf_gas_mix_eos : therm_perf_gas_mix_eos.d $(LIBLUA) $(NM_FILES)
	dmd -main -unittest therm_perf_gas_mix_eos.d ../gas_model.d \
		../physical_constants.d evt_eos.d cea_thermo_curves.d \
		$(UTIL_FILES) $(NM_FILES) -I../../ -I$(UTIL_DIR) $(DLINKFLAGS)
