# makefile for the gas module

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd
PLATFORM ?= linux
INSTALL_DIR ?= $(HOME)/dgdinst
BUILD_DIR := ../../build

DEMO_PROGRAMS := gas_model_demo ideal_gas_demo co2_model_demo co2_sw_model_demo \
	sf6_model_demo luagas_model_demo univariate_lut_demo build_tree 
TEST_PROGRAMS := cea_therm_cond_test cea_thermo_curves_test cea_viscosity_test \
	gas_model_test ideal_gas_test powers_aslam_gas_test \
	perf_gas_mix_eos_test sutherland_therm_cond_test sutherland_viscosity_test \
	therm_perf_gas_test therm_perf_gas_mix_eos_test very_viscous_air_test \
	co2gas_sw_test uniform_lut_test adaptive_lut_CEA_test ideal_air_proxy_test

LUA := ../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib
ifeq ($(DMD), dmd)
    # DFLAGS := -w
    DFLAGS := -w -O -release -inline -boundscheck=off
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -O -release -inline -boundscheck=off
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), gdc)
    DFLAGS := -Wall -O2 -frelease -fno-debug
    OF := -o
    DLINKFLAGS := -L$(LIBLUAPATH) $(LIBLUA) -ldl
endif

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

include gas_files.mk

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm *.o
	- rm ideal_air_fortran.mod
	- rm ${DEMO_PROGRAMS}
	- rm ${TEST_PROGRAMS}
	- rm gas.so
	- rm gas-calc
	- rm species-database/species-database.lua
	- rm LOGFILE*TEST
	- cd $(LUA); make clean

test: $(TEST_PROGRAMS) gas-calc
	tclsh gas-package-test.tcl

build-prep-gas: prep_gas.lua species-database.lua
	- mkdir -p $(BUILD_DIR)/bin
	cp prep_gas.lua $(BUILD_DIR)/bin/prep-gas; chmod +x $(BUILD_DIR)/bin/prep-gas
	- mkdir -p $(BUILD_DIR)/data
	cp species-database/species-database.lua $(BUILD_DIR)/data/
	cp species-database/species-list.txt $(BUILD_DIR)/data/

gas-calc: luagas_model.d $(LIBLUA) $(GAS_FILES) $(GAS_LUA_FILES) $(NM_FILES) \
	$(UTIL_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES)
	dmd -ofgas-calc -version=gas_calc \
		$(GAS_FILES) $(GAS_LUA_FILES) $(NM_FILES) $(UTIL_FILES) $(KINETICS_FILES) \
		$(KINETICS_LUA_FILES) $(DLINKFLAGS)

species-data-converter: species_data_converter.lua
	- mkdir -p $(BUILD_DIR)/bin
	cp species_data_converter.lua $(BUILD_DIR)/bin/species-data-converter; \
		chmod +x $(BUILD_DIR)/bin/species-data-converter

install: gas-calc build-prep-gas species-data-converter
	- mkdir -p $(INSTALL_DIR)
	cp gas-calc $(BUILD_DIR)/bin/
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

$(LIBLUA):
	cd $(LUA); make linux local

species-database.lua :
	cd species-database; make species-database.lua

gas_model_demo: gas_model_demo.d $(LIBLUA) 
	$(DMD) -g $(DFLAGS) gas_model_demo.d $(OF)gas_model_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS) 

ideal_gas_demo: ideal_gas_demo.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) ideal_gas_demo.d $(OF)ideal_gas_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

co2_model_demo: co2_model_demo.d  $(LIBLUA) $(UTIL_FILES) $(NM_FILES) 
	$(DMD) -g $(DFLAGS) co2_model_demo.d $(OF)co2_model_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

co2_sw_model_demo: co2_sw_model_demo.d  $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) co2_sw_model_demo.d $(OF)co2_sw_model_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

sf6_model_demo: sf6_model_demo.d  $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) sf6_model_demo.d $(OF)sf6_model_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

luagas_model_demo: luagas_model_demo.d $(GAS_LUA_FILES) $(LIBLUA) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES)
	$(DMD) -g $(DFLAGS) luagas_model_demo.d $(OF)luagas_model_demo \
		$(GAS_FILES) $(GAS_LUA_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS) 


# Jonathan's programs that were moved out of nm.

build_tree: build_tree.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) build_tree.d $(OF)build_tree \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

univariate_lut_demo: univariate_lut_demo.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -debug univariate_lut_demo.d $(OF)univariate_lut_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

#ideal_gas_test : ideal_gas.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
#	$(DMD) -main -unittest $(DFLAGS) $(OF)ideal_gas_test \
#		$(GAS_FILES) $(UTIL_FILES) $(NM_FILES) \
#		$(DLINKFLAGS)

#therm_perf_gas_test : therm_perf_gas.d libgas.a $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
#	$(DMD) -main -unittest $(DFLAGS) therm_perf_gas.d $(OF)therm_perf_gas_test \
#		libgas.a $(UTIL_FILES) $(NM_FILES) \
#		$(DLINKFLAGS)

#very_viscous_air_test : very_viscous_air.d libgas.a $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
#	$(DMD) -main -unittest $(DFLAGS) very_viscous_air.d $(OF)very_viscous_air_test \
#		libgas.a $(UTIL_FILES) $(NM_FILES) \
#		$(DLINKFLAGS)

cea_therm_cond_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

cea_viscosity_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

cea_thermo_curves_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

gas_model_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

ideal_gas_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

perf_gas_mix_eos_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

sutherland_therm_cond_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

sutherland_viscosity_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

therm_perf_gas_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

therm_perf_gas_mix_eos_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

very_viscous_air_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

co2gas_sw_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

uniform_lut_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

adaptive_lut_CEA_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

lut_comp_demo: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA) lut_comp_demo.d
	dmd -of$@ -debug -gc -w lut_comp_demo.d \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

ideal_air_proxy_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

powers_aslam_gas_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	dmd -of$@ -debug -gc -w -version=$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

#--------------------------------------------------------------------------
# Fortran bits

ideal_air_fortran.o: ideal_air_fortran.f
	gfortran -c -ffree-form ideal_air_fortran.f
