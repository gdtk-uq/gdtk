# makefile for the gas module

DMD ?= ldc2

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug
WITH_FPE ?= 1

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR := ../../build

DEMO_PROGRAMS := gas_model_demo ideal_gas_demo luagas_model_demo

TEST_PROGRAMS := chemkin_therm_cond_test chemkin_viscosity_test cea_therm_cond_test \
	cea_thermo_curves_test cea_viscosity_test composite_gas_test \
	init_gas_model_test ideal_gas_test ideal_helium_test \
	cubic_gas_test cea_gas_test ideal_gas_ab_test \
	ideal_dissociating_gas_test fuel_air_mix_test equilibrium_gas_test \
	perf_gas_mix_eos_test sutherland_therm_cond_test sutherland_viscosity_test \
	wilke_mixing_therm_cond_test wilke_mixing_viscosity_test \
	therm_perf_gas_test therm_perf_gas_complex_test therm_perf_gas_equil_test \
	two_temperature_gas_test two_temperature_trans_props_test \
	therm_perf_gas_mix_eos_test very_viscous_air_test \
	uniform_lut_test uniform_lut_plus_ideal_test \
	adaptive_lut_CEA_test \
        two_temperature_reacting_argon_test two_temperature_argon_plus_ideal_test \
	init_gas_model_complex_test ideal_gas_complex_test cea_thermo_curves_complex_test \
	vib_specific_nitrogen_test vib_specific_co_test \
	vib_specific_co_mixture_test \
	two_temperature_gasgiant_test multi_temperature_gas_test

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib

CEQ_DIR := ../extern/ceq/source
LIBCEQ := $(CEQ_DIR)/libceq.a
include $(CEQ_DIR)/ceq_files.mk

NTYPES_DIR := ../ntypes
NTYPES_FILES := $(NTYPES_DIR)/complex.d

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
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
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
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
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
        DFLAGS := -w -g -O -release -enable-inlining -boundscheck=off -d-version=flavour_fast -ffast-math
        ifeq ($(PLATFORM), linux)
	    # Full link-time optimization does not play nicely on macOS
	    DFLAGS += -flto=full
        endif
    endif
    PIC := --relocation-model=pic
    DVERSION := -d-version=
    OF := -of=
    DLINKFLAGS := -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif

ifeq ($(WITH_FPE),1)
    DFLAGS += $(DVERSION)enable_fp_exceptions
endif

# DIP1008 allows throwing of exceptions in @nogc code.
# See notes in src/eilmer/makefile.
DFLAGS += -dip1008 -preview=in

# ----------------------------------------------------------------------
# Here begins the list of targets, starting with the top-level actions.
# ----------------------------------------------------------------------

CIF_DIR := species-database/collision-integrals
COLLISION_INTEGRAL_FILES := $(CIF_DIR)/gupta_etal_1990_CI_data.lua \
	$(CIF_DIR)/wright_etal_CI_data.lua \
	$(CIF_DIR)/palmer_etal_CI_data.lua

.PHONY: build-prep-gas
build-prep-gas: $(BUILD_DIR)/bin/prep-gas

# prep-gas and its database is required by the Eilmer flow solver.
$(BUILD_DIR)/bin/prep-gas: prep_gas.lua $(BUILD_DIR)/data/species-database.lua \
		species-database/species-list.txt species_data_converter.lua $(COLLISION_INTEGRAL_FILES)
	@mkdir -p $(BUILD_DIR)/bin
	cp prep_gas.lua $(BUILD_DIR)/bin/prep-gas; chmod +x $(BUILD_DIR)/bin/prep-gas
	cp species_data_converter.lua $(BUILD_DIR)/bin/species-data-converter; \
		chmod +x $(BUILD_DIR)/bin/species-data-converter
	@mkdir -p $(BUILD_DIR)/data
	cp species-database/species-list.txt $(BUILD_DIR)/data/
	cp $(COLLISION_INTEGRAL_FILES) $(BUILD_DIR)/data/

$(BUILD_DIR)/data/species-database.lua:
	@mkdir $(BUILD_DIR)/data
	$(MAKE) SPECIES_DATABASE_DIR=$(BUILD_DIR)/data -C species-database
