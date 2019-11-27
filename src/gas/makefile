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
	cubic_gas_test cea_gas_test powers_aslam_gas_test \
	ideal_dissociating_gas_test fuel_air_mix_test equilibrium_gas_test \
	perf_gas_mix_eos_test sutherland_therm_cond_test sutherland_viscosity_test \
	wilke_mixing_therm_cond_test wilke_mixing_viscosity_test \
	therm_perf_gas_test therm_perf_gas_mix_eos_test very_viscous_air_test \
	co2gas_sw_test uniform_lut_test uniform_lut_plus_ideal_test \
	adaptive_lut_CEA_test ideal_air_proxy_test \
        two_temperature_reacting_argon_test \
	steam_test \
	init_gas_model_complex_test ideal_gas_complex_test cea_thermo_curves_complex_test \
	pseudo_species_gas_test vib_specific_nitrogen_test \
	electronically_specific_gas_test \
	two_temperature_gasgiant_test

LUA := ../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib
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
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
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
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
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
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), gdc)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -Wall -Og -g -version=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        DFLAGS := -Wall -O2 -frelease -fno-debug -version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -Wall -O2 -frelease -fno-debug -version=flavour_fast
    endif
    PIC := -fPIC
    DVERSION := -version=
    OF := -o
    DLINKFLAGS := -L$(LIBLUAPATH) $(LIBLUA) -ldl
endif

ifeq ($(WITH_DVODE),1)
    DFLAGS += -version=with_dvode
    DLINKFLAGS += -L-lgfortran
endif

# DIP1008 allows throwing of exceptions in @nogc code.
DFLAGS += -dip1008

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

include gas_files.mk

# ----------------------------------------------------------------------
# Here begins the list of targets, starting with the top-level actions.
# ----------------------------------------------------------------------

build: gas-calc prep_gas.lua species-database/species-database.lua \
		species-database/species-list.txt species_data_converter.lua \
		libgas.so 
	- mkdir -p $(BUILD_DIR)/bin
	cp gas-calc $(BUILD_DIR)/bin/
	cp prep_gas.lua $(BUILD_DIR)/bin/prep-gas; chmod +x $(BUILD_DIR)/bin/prep-gas
	cp species_data_converter.lua $(BUILD_DIR)/bin/species-data-converter; \
		chmod +x $(BUILD_DIR)/bin/species-data-converter
	- mkdir -p $(BUILD_DIR)/data
	cp species-database/species-database.lua $(BUILD_DIR)/data/
	cp species-database/species-list.txt $(BUILD_DIR)/data/
	- mkdir -p $(BUILD_DIR)/lib
	cp libgas.so $(BUILD_DIR)/lib/
	cp -r ../lib/* $(BUILD_DIR)/lib/

build-prep-gas: prep_gas.lua species-database/species-database.lua \
		species-database/species-list.txt species_data_converter.lua
	- mkdir -p $(BUILD_DIR)/bin
	cp prep_gas.lua $(BUILD_DIR)/bin/prep-gas; chmod +x $(BUILD_DIR)/bin/prep-gas
	cp species_data_converter.lua $(BUILD_DIR)/bin/species-data-converter; \
		chmod +x $(BUILD_DIR)/bin/species-data-converter
	- mkdir -p $(BUILD_DIR)/data
	cp species-database/species-database.lua $(BUILD_DIR)/data/
	cp species-database/species-list.txt $(BUILD_DIR)/data/

install: build
	- mkdir -p $(INSTALL_DIR)
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

test: $(TEST_PROGRAMS) gas-calc
	tclsh gas-package-test.tcl

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm *.o *.so
	- rm *.mod
	- rm -r $(BUILD_DIR)/*
	- rm -r __pycache__
	- rm ${DEMO_PROGRAMS}
	- rm ${TEST_PROGRAMS}
	- rm gas.so
	- rm gas-calc
	- rm species-database/species-database.lua
	- rm LOGFILE*TEST
	- rm thermo.inp thermo.lib thermo.out
	- rm trans.inp trans.lib trans.out
	- rm tmp.inp tmp.out tmp.plt
	- rm ./ideal-air-gas-model.lua ./cea-lut-air-version-test.lua
	- rm ./uniform-lut-plus-ideal-air-gas-model.lua
	- cd $(LUA); make clean
	- rm $(KINETICS_DIR)/pseudosp_solve_ode.o
	- rm $(KINETICS_DIR)/pseudosp_rhs_ode.o
	- rm $(KINETICS_DIR)/dvode.o
	- rm $(KINETICS_DIR)/updatepseudospecieschem.mod

# Specific targets for individual artefacts.

gas-calc: luagas_model.d $(LIBLUA) $(GAS_FILES) $(GAS_LUA_FILES) $(NM_FILES) \
	$(UTIL_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES)
	$(DMD) $(OF)gas-calc $(DVERSION)gas_calc $(DFLAGS) \
		$(GAS_FILES) $(GAS_LUA_FILES) $(NM_FILES) $(UTIL_FILES) $(KINETICS_FILES) \
		$(KINETICS_LUA_FILES) $(DLINKFLAGS)

libgas.so: gas_cwrap.d $(LIBLUA) $(GAS_FILES) $(GAS_LUA_FILES) $(NM_FILES) \
	$(UTIL_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES)
	$(DMD) -c gas_cwrap.d $(PIC) $(DFLAGS) \
		$(GAS_FILES) $(NM_FILES) $(UTIL_FILES) $(KINETICS_FILES)
	$(DMD) $(OF)libgas.so gas_cwrap.o $(DFLAGS) \
		-shared -defaultlib=libphobos2.so \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		-L-rpath=$(INSTALL_DIR)/lib $(DLINKFLAGS)

$(LIBLUA):
	cd $(LUA); make linux local

species-database/species-database.lua :
	cd species-database; make species-database.lua

# Some demo targets. These are left over from the early days of development.

gas_model_demo: gas_model_demo.d $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) -g $(DFLAGS) gas_model_demo.d $(OF)gas_model_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

ideal_gas_demo: ideal_gas_demo.d $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) ideal_gas_demo.d $(OF)ideal_gas_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

co2_model_demo: co2_model_demo.d $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) co2_model_demo.d $(OF)co2_model_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

co2_sw_model_demo: co2_sw_model_demo.d $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) co2_sw_model_demo.d $(OF)co2_sw_model_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

sf6_model_demo: sf6_model_demo.d $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g $(DFLAGS) sf6_model_demo.d $(OF)sf6_model_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

luagas_model_demo: luagas_model_demo.d $(GAS_FILES) $(GAS_LUA_FILES) $(LIBLUA) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES)
	$(DMD) -g $(DFLAGS) luagas_model_demo.d $(OF)luagas_model_demo \
		$(GAS_FILES) $(GAS_LUA_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

# Jonathan's programs that were moved out of nm.

build_tree: build_tree.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) build_tree.d -dip1008 $(OF)build_tree \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

univariate_lut_demo: univariate_lut_demo.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -debug -dip1008 univariate_lut_demo.d $(OF)univariate_lut_demo \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(DLINKFLAGS)

# Individual test targets

cea_therm_cond_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

chemkin_therm_cond_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

chemkin_viscosity_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

cea_viscosity_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

cea_thermo_curves_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

cea_thermo_curves_complex_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)cea_thermo_curves_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)


init_gas_model_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

init_gas_model_complex_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)init_gas_model_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

ideal_gas_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

ideal_gas_complex_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)ideal_gas_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

ideal_helium_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

cubic_gas_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

cea_gas_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

perf_gas_mix_eos_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

sutherland_therm_cond_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

sutherland_viscosity_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

therm_perf_gas_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

therm_perf_gas_mix_eos_test: $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

very_viscous_air_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

wilke_mixing_therm_cond_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

wilke_mixing_viscosity_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

co2gas_sw_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

uniform_lut_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

uniform_lut_plus_ideal_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

adaptive_lut_CEA_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

lut_comp_demo: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA) lut_comp_demo.d
	$(DMD) $(OF)$@ $(DFLAGS) lut_comp_demo.d \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

ideal_air_proxy_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

powers_aslam_gas_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

vib_specific_nitrogen_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

two_temperature_reacting_argon_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

ideal_dissociating_gas_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

fuel_air_mix_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

equilibrium_gas_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

steam_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

pseudo_species_gas_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

electronically_specific_gas_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

two_temperature_gasgiant_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

gasgiant_transport_properties_test: $(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

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
