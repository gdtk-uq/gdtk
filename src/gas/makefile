# makefile for the gas module

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd
INSTALL_DIR ?= $(HOME)/e3bin

DEMO_PROGRAMS := gas_model_demo ideal_gas_demo CO2_model_demo
UNITTEST_FILES := therm_perf_gas_test \
	ideal_gas_test

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua.d \
	$(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

NM_DIR := ../nm
NM_FILES := $(NM_DIR)/ridder.d

LUA := ../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib

GAS_DIR := .
GAS_MODEL_FILES := \
	$(GAS_DIR)/package.d \
	$(GAS_DIR)/gas_model.d \
	$(GAS_DIR)/gas_model_util.d \
	$(GAS_DIR)/ideal_gas.d \
	$(GAS_DIR)/physical_constants.d \
	$(GAS_DIR)/therm_perf_gas.d \
	$(GAS_DIR)/very_viscous_air.d\
	$(GAS_DIR)/CO2Gas.d

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

ifeq ($(DMD), dmd)
    # DFLAGS := -w
    DFLAGS := -w -O -release -inline -boundscheck=off \
	-I.. -I../util -I../extern/LuaD -I../nm
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
    LIBNAME := libgas
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -O -release -inline -boundscheck=off \
	-I.. -I../util -I../nm
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
    LIBNAME := ./gas
endif
ifeq ($(DMD), gdc)
    DFLAGS := -Wall -O2 -frelease -fno-debug \
	-I.. -I../util -I../extern/LuaD -I../nm
    OF := -o
    DLINKFLAGS := -L$(LIBLUAPATH) $(LIBLUA) -ldl
endif

libgas.a : $(GAS_FILES) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -lib $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(OF)$(LIBNAME) $(DFLAGS)

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm ${DEMO_PROGRAMS} *.o *.a
	- cd $(LUA); make clean

install-prep-gas: prep_gas.lua
	cp prep_gas.lua $(INSTALL_DIR)/prep-gas; chmod +x $(INSTALL_DIR)/prep-gas

$(LIBLUA):
	cd $(LUA); make linux local

gas_model_demo: gas_model_demo.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) gas_model_demo.d $(OF)gas_model_demo \
		$(GAS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

ideal_gas_demo: ideal_gas_demo.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) ideal_gas_demo.d $(OF)ideal_gas_demo \
		$(GAS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)
		
CO2_model_demo: CO2_model_demo.d  $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) CO2_model_demo.d $(OF)CO2_model_demo \
		$(GAS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

ideal_gas_test : ideal_gas.d libgas.a $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -main -unittest $(DFLAGS) ideal_gas.d $(OF)ideal_gas_test \
		libgas.a $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

therm_perf_gas_test : therm_perf_gas.d libgas.a $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -main -unittest $(DFLAGS) therm_perf_gas.d $(OF)therm_perf_gas_test \
		libgas.a $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

very_viscous_air_test : very_viscous_air.d libgas.a $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -main -unittest $(DFLAGS) very_viscous_air.d $(OF)very_viscous_air_test \
		libgas.a $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)
