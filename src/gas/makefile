# makefile for the gas module

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug

PLATFORM ?= linux
WITH_DVODE ?= 0
INSTALL_DIR ?= $(HOME)/dgdinst
BUILD_DIR := ../../build

DEMO_PROGRAMS := gas_model_demo ideal_gas_demo co2_model_demo co2_sw_model_demo \
	sf6_model_demo luagas_model_demo univariate_lut_demo build_tree
TEST_PROGRAMS := chemkin_therm_cond_test chemkin_viscosity_test cea_therm_cond_test \
	cea_thermo_curves_test cea_viscosity_test \
	init_gas_model_test ideal_gas_test ideal_helium_test \
	cubic_gas_test cea_gas_test ideal_gas_ab_test \
	ideal_dissociating_gas_test fuel_air_mix_test equilibrium_gas_test \
	perf_gas_mix_eos_test sutherland_therm_cond_test sutherland_viscosity_test \
	wilke_mixing_therm_cond_test wilke_mixing_viscosity_test \
	therm_perf_gas_test therm_perf_gas_complex_test therm_perf_gas_equil_test \
	therm_perf_gas_mix_eos_test very_viscous_air_test \
	co2gas_sw_test uniform_lut_test uniform_lut_plus_ideal_test \
	adaptive_lut_CEA_test ideal_air_proxy_test \
        two_temperature_reacting_argon_test two_temperature_argon_plus_ideal_test \
	steam_test \
	init_gas_model_complex_test ideal_gas_complex_test cea_thermo_curves_complex_test \
	pseudo_species_gas_test vib_specific_nitrogen_test \
	electronically_specific_gas_test \
	two_temperature_gasgiant_test

LUA := ../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib

CEQ_DIR := ../extern/ceq/source
LIBCEQ := $(CEQ_DIR)/libceq.a
include $(CEQ_DIR)/ceq_files.mk

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

GASDYN_DIR := ../gasdyn
include $(GASDYN_DIR)/gasdyn_files.mk

GEOM_DIR := ../geom
include $(GEOM_DIR)/geom_files.mk

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

include gas_files.mk

ifeq ($(DMD), dmd)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -w -g -debug -version=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        DFLAGS := -profile -w -g -O -release -boundscheck=off -version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O -release -boundscheck=off -version=flavour_fast
    endif
    PIC := -fPIC
    DVERSION := -version=
    OF := -of
    DLINKFLAGS := -L-ldl
endif
ifeq ($(DMD), ldmd2)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -w -g -debug -version=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        DFLAGS := -profile -w -g -O -release -inline -boundscheck=off -version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O -release -inline -boundscheck=off -version=flavour_fast
    endif
    PIC := -fPIC
    DVERSION := -version=
    OF := -of
    DLINKFLAGS := -L-ldl
endif
ifeq ($(DMD), ldc2)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -w -g -d-debug -d-version=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        # -fprofile-generate will result in profraw files being written
        # that may be viewed, showing the top 10 functions with internal block counts
        # llvm-profdata show -text -topn=10 <profraw-file>
        DFLAGS := -fprofile-generate -g -w -O -release -enable-inlining -boundscheck=off -d-version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O -release -enable-inlining -boundscheck=off -d-version=flavour_fast -ffast-math -flto=full
    endif
    PIC := --relocation-model=pic
    DVERSION := -d-version=
    OF := -of=
    DLINKFLAGS := -L-ldl
endif
ifeq ($(DMD), gdc-8)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -Wall -Og -g -fversion=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        DFLAGS := -Wall -O2 -frelease -fno-debug -fversion=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -Wall -O2 -frelease -fno-debug -fversion=flavour_fast
    endif
    PIC := -fPIC
    DVERSION := -fversion=
    OF := -o
    DLINKFLAGS := -ldl
endif

ifeq ($(WITH_DVODE),1)
    DFLAGS += -version=with_dvode
    DLINKFLAGS += -L-lgfortran
endif

# DIP1008 allows throwing of exceptions in @nogc code.
# See notes in src/eilmer/makefile.
DFLAGS += -dip1008

# ----------------------------------------------------------------------
# Here begins the list of targets, starting with the top-level actions.
# ----------------------------------------------------------------------

# prep-gas and its database is required by the Eilmer flow solver.
build-prep-gas: prep_gas.lua species-database/species-database.lua \
		species-database/species-list.txt species_data_converter.lua
	- mkdir -p $(BUILD_DIR)/bin
	cp prep_gas.lua $(BUILD_DIR)/bin/prep-gas; chmod +x $(BUILD_DIR)/bin/prep-gas
	cp species_data_converter.lua $(BUILD_DIR)/bin/species-data-converter; \
		chmod +x $(BUILD_DIR)/bin/species-data-converter
	- mkdir -p $(BUILD_DIR)/data
	cp species-database/species-database.lua $(BUILD_DIR)/data/
	cp species-database/species-list.txt $(BUILD_DIR)/data/

# The following are not required by the Eilmer flow solver but
# they are useful for building stand-alone flow-analysis programs.
build-libgas: libgas.so gas-calc build-uniform-lut.py $(LIBCEQ)
	- mkdir -p $(BUILD_DIR)/bin
	cp gas-calc $(BUILD_DIR)/bin/
	cp build-uniform-lut.py $(BUILD_DIR)/bin/build-uniform-lut; chmod +x $(BUILD_DIR)/bin/build-uniform-lut
	cp ../../examples/estcn/estcn.py $(BUILD_DIR)/bin/estcn; chmod +x $(BUILD_DIR)/bin/estcn
	- mkdir -p $(BUILD_DIR)/lib
	cp libgas.so $(BUILD_DIR)/lib/
	cp $(CEQ_DIR)/libceq.so $(BUILD_DIR)/lib/
	cp -r ../lib/* $(BUILD_DIR)/lib/

install: build-prep-gas build-libgas
	- mkdir -p $(INSTALL_DIR)
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

test: $(TEST_PROGRAMS) gas-calc
	tclsh gas-package-test.tcl

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm -f *.o *.so
	- rm -f *.mod
	- rm -rf $(BUILD_DIR)/*
	- rm -rf __pycache__
	- rm -f ${DEMO_PROGRAMS}
	- rm -f ${TEST_PROGRAMS}
	- rm -f gas.so
	- rm -f gas-calc
	- rm -f species-database/species-database.lua
	- rm -f LOGFILE*TEST
	- rm -f thermo.inp thermo.lib thermo.out
	- rm -f trans.inp trans.lib trans.out
	- rm -f tmp.inp tmp.out tmp.plt
	- rm -f air-eq-1.lua air-eq.inp air-eq.lua
	- rm -f ./ideal-air-gas-model.lua ./cea-lut-air-version-test.lua
	- rm -f ./uniform-lut-plus-ideal-air-gas-model.lua
	- rm -f ./two-temperature-reacting-argon-model.lua
	- rm -f ./two-temperature-argon-plus-ideal-air-gas-model.lua
	- rm -f ./diffusion/*.o
	- rm -f ./thermo/*.o
	- rm -f $(GASDYN_DIR)/*.o
	- cd $(LUA); make clean
	- cd $(CEQ_DIR); make clean
	- cd $(UTIL_DIR); make clean
	- cd $(GEOM_DIR); make clean
	- cd $(GZIP_DIR); make clean
	- cd $(NM_DIR); make clean
	- cd $(KINETICS_DIR); make clean

# Specific targets for individual artefacts.

gas-calc: luagas_model.d $(LIBLUA) $(GAS_FILES) $(GAS_LUA_FILES) \
	$(NM_FILES) $(NM_LUA_FILES) $(UTIL_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(GEOM_FILES) $(GZIP_FILES) \
	$(CEQ_SRC_FILES) $(LIBCEQ)
	$(DMD) -c $(OF)gas_calc.o $(DVERSION)gas_calc $(DFLAGS) \
		$(GAS_FILES) $(GAS_LUA_FILES) \
		$(NM_FILES) $(NM_LUA_FILES) $(UTIL_FILES) \
		$(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(CEQ_SRC_FILES) $(GEOM_FILES) $(GZIP_FILES)
	$(DMD) $(OF)gas-calc $(DFLAGS) gas_calc.o \
		ideal_air_fortran.o $(LIBCEQ) $(LIBLUA) \
		$(DLINKFLAGS)

libgas.so: gas_cwrap.d $(LIBLUA) $(GAS_FILES) $(KINETICS_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(GASDYN_FILES) \
	$(CEQ_SRC_FILES) $(LIBCEQ) $(GEOM_FILES) $(GZIP_FILES)
	$(DMD) -c gas_cwrap.d $(PIC) $(DFLAGS) -op \
		$(GAS_FILES) $(KINETICS_FILES) $(GASDYN_FILES) $(NM_FILES) \
		$(UTIL_FILES) $(CEQ_SRC_FILES) $(GEOM_FILES) $(GZIP_FILES)
	$(DMD) $(OF)libgas.so gas_cwrap.o $(PIC) $(DFLAGS) -op -shared \
		$(GAS_FILES) $(KINETICS_FILES) $(GASDYN_FILES) $(NM_FILES) \
		$(UTIL_FILES) $(CEQ_SRC_FILES) $(GEOM_FILES) $(GZIP_FILES) \
		$(LIBCEQ) $(LIBLUA) \
		-L-L=$(INSTALL_DIR)/lib $(DLINKFLAGS)

$(LIBLUA):
	cd $(LUA); make $(PLATFORM) local

$(LIBCEQ):
	cd $(CEQ_DIR); make

species-database/species-database.lua :
	cd species-database; make species-database.lua

# Some demo targets. These are left over from the early days of development.

gas_model_demo: gas_model_demo.d $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) -g $(DFLAGS) gas_model_demo.d $(OF)gas_model_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA) \
		$(DLINKFLAGS)

ideal_gas_demo: ideal_gas_demo.d $(GAS_FILES) $(CEQ_SRC_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) ideal_gas_demo.d $(OF)ideal_gas_demo \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA) \
		$(DLINKFLAGS)

co2_model_demo: co2_model_demo.d $(GAS_FILES) $(CEQ_SRC_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) co2_model_demo.d $(OF)co2_model_demo \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA) \
		$(DLINKFLAGS)

co2_sw_model_demo: co2_sw_model_demo.d $(GAS_FILES) $(CEQ_SRC_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) co2_sw_model_demo.d $(OF)co2_sw_model_demo \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA) \
		$(DLINKFLAGS)

sf6_model_demo: sf6_model_demo.d $(GAS_FILES) $(CEQ_SRC_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) sf6_model_demo.d $(OF)sf6_model_demo \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA) \
		$(DLINKFLAGS)

luagas_model_demo: luagas_model_demo.d $(GAS_FILES) $(CEQ_SRC_FILES) $(GAS_LUA_FILES) $(LIBLUA) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES)
	$(DMD) -g $(DFLAGS) luagas_model_demo.d $(OF)luagas_model_demo \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GAS_LUA_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(LIBLUA) \
		$(DLINKFLAGS)

# Jonathan's programs that were moved out of nm.

build_tree: build_tree.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) build_tree.d -dip1008 $(OF)build_tree \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA) \
		$(DLINKFLAGS)

univariate_lut_demo: univariate_lut_demo.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -debug -dip1008 univariate_lut_demo.d $(OF)univariate_lut_demo \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA) \
		$(DLINKFLAGS)

# Individual test targets

cea_therm_cond_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

chemkin_therm_cond_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

chemkin_viscosity_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

cea_viscosity_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

cea_thermo_curves_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

cea_thermo_curves_complex_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)cea_thermo_curves_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)


init_gas_model_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

init_gas_model_complex_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)init_gas_model_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

ideal_gas_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

ideal_gas_complex_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)ideal_gas_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

ideal_helium_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

cubic_gas_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

cea_gas_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

perf_gas_mix_eos_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

sutherland_therm_cond_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

sutherland_viscosity_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

therm_perf_gas_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

therm_perf_gas_complex_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)therm_perf_gas_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

therm_perf_gas_equil_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

therm_perf_gas_mix_eos_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

very_viscous_air_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

wilke_mixing_therm_cond_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA) $(LIBCEQ) $(DLINKFLAGS)

wilke_mixing_viscosity_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

co2gas_sw_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

uniform_lut_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

uniform_lut_plus_ideal_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

adaptive_lut_CEA_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

lut_comp_demo: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) lut_comp_demo.d
	$(DMD) $(OF)$@ $(DFLAGS) lut_comp_demo.d \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

ideal_air_proxy_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

ideal_gas_ab_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

vib_specific_nitrogen_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

two_temperature_reacting_argon_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

two_temperature_argon_plus_ideal_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

ideal_dissociating_gas_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

fuel_air_mix_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

equilibrium_gas_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

steam_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

pseudo_species_gas_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

electronically_specific_gas_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

two_temperature_gasgiant_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

gasgiant_transport_properties_test: $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

# Fortran bits

ideal_air_fortran.o: ideal_air_fortran.f
	gfortran -c -ffree-form -fPIC ideal_air_fortran.f

FFLAGS_fort=-freal-4-real-8 -fdefault-real-8  -fdefault-double-8
FFLAGS_DEBUG_fort=-O0 -g -ggdb -finit-real=snan -ffpe-trap=zero,overflow -fcheck=all -std=legacy -Wall -fbacktrace

$(KINETICS_DIR)/pseudosp_solve_ode.o: $(KINETICS_DIR)/pseudosp_solve_ode.f90
	cd $(KINETICS_DIR); gfortran -c $(FFLAGS_fort) $(FFLAGS_DEBUG_fort) pseudosp_solve_ode.f90

$(KINETICS_DIR)/pseudosp_rhs_ode.o: $(KINETICS_DIR)/pseudosp_rhs_ode.f90
	cd $(KINETICS_DIR); gfortran -c $(FFLAGS_fort) $(FFLAGS_DEBUG_fort) pseudosp_rhs_ode.f90

$(KINETICS_DIR)/dvode.o: $(KINETICS_DIR)/dvode.f
	cd $(KINETICS_DIR); gfortran -c $(FFLAGS_fort) dvode.f
