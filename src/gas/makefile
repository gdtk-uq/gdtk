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
PY_PROGRAMS := build-uniform-lut

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

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a
include $(EQC_DIR)/eqc_files.mk

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

GDTK_PY_DIR := ../gdtk-py

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
        DFLAGS := -w -g -link-defaultlib-debug -d-debug -d-version=flavour_debug
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

# prep-gas and its database is required by the Eilmer flow solver.
build-prep-gas: prep_gas.lua species-database/species-database.lua \
		species-database/species-list.txt species_data_converter.lua $(COLLISION_INTEGRAL_FILES)
	- mkdir -p $(BUILD_DIR)/bin
	cp prep_gas.lua $(BUILD_DIR)/bin/prep-gas; chmod +x $(BUILD_DIR)/bin/prep-gas
	cp species_data_converter.lua $(BUILD_DIR)/bin/species-data-converter; \
		chmod +x $(BUILD_DIR)/bin/species-data-converter
	- mkdir -p $(BUILD_DIR)/data
	cp species-database/species-database.lua $(BUILD_DIR)/data/
	cp species-database/species-list.txt $(BUILD_DIR)/data/
	cp $(COLLISION_INTEGRAL_FILES) $(BUILD_DIR)/data/

# The following are not required by the Eilmer flow solver but
# they are useful for building stand-alone flow-analysis programs.
build-libgas: libgas.so gas-calc build-uniform-lut $(LIBEQC)
	- mkdir -p $(BUILD_DIR)/bin
	cp gas-calc $(BUILD_DIR)/bin/
	cp ../../examples/estcn/estcn.py $(BUILD_DIR)/bin/estcn; chmod +x $(BUILD_DIR)/bin/estcn
	- mkdir -p $(BUILD_DIR)/lib
	cp -r ../lib/* $(BUILD_DIR)/lib/
	make -C $(GDTK_PY_DIR) build
	cp -r $(GDTK_PY_DIR)/gdtk $(BUILD_DIR)/lib/

install: build-prep-gas build-libgas
	- mkdir -p $(INSTALL_DIR)
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

test: $(TEST_PROGRAMS) gas-calc
	tclsh gas-package-test.tcl

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm -f *.o *.so *.a *.obj *.lib
	- rm -f *.mod
	- rm -rf $(BUILD_DIR)/*
	- rm -rf __pycache__/ *.egg-info/
	- rm -f ${DEMO_PROGRAMS}
	- rm -f ${TEST_PROGRAMS}
	- rm -f gas-calc
	- rm -f species-database/species-database.lua
	- rm -f LOGFILE*TEST
	- rm -f thermo.inp thermo.lib thermo.out
	- rm -f trans.inp trans.lib trans.out
	- rm -f tmp.inp tmp.out tmp.plt
	- rm -f air-eq-1.lua air-eq.inp air-eq.lua
	- rm -f air-5sp-1T.inp air-5sp-1T.lua air-5sp-eq.lua
	- rm -f ./ideal-air-gas-model.lua ./cea-lut-air-version-test.lua
	- rm -f ./uniform-lut-plus-ideal-air-gas-model.lua
	- rm -f ./two-temperature-reacting-argon-model.lua
	- rm -f ./two-temperature-argon-plus-ideal-air-gas-model.lua
	- rm -f ./diffusion/*.o ./diffusion/*.obj
	- rm -f ./thermo/*.o ./thermo/*.obj
	- rm -f $(GASDYN_DIR)/*.o $(GASDYN_DIR)/*.obj
	- cd $(LUA_DIR); make clean
	- cd $(EQC_DIR); make clean
	- cd $(UTIL_DIR); make clean
	- cd $(GEOM_DIR); make clean
	- cd $(GZIP_DIR); make clean
	- cd $(NM_DIR); make clean
	- cd $(KINETICS_DIR); make clean

# Specific targets for individual artefacts.

gas-calc: gas_calc.d $(LIBLUA) $(GAS_FILES) $(GAS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(NM_LUA_FILES) $(UTIL_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(GEOM_FILES) $(GZIP_FILES) \
	$(EQC_SRC_FILES) $(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)gas-calc $(DFLAGS) gas_calc.d \
		$(GAS_FILES) $(GAS_LUA_FILES) \
		$(NM_FILES) $(NTYPES_FILES) $(NM_LUA_FILES) $(UTIL_FILES) \
		$(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(EQC_SRC_FILES) $(GEOM_FILES) $(GZIP_FILES) \
		$(LIBEQC) $(LIBLUA) \
		$(DLINKFLAGS)

libgas.so: gas_cwrap.d $(LIBLUA) $(GAS_FILES) $(GAS_LUA_FILES) $(KINETICS_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(GASDYN_FILES) \
	$(EQC_SRC_FILES) $(GEOM_FILES) $(GZIP_FILES) $(KINETICS_LUA_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) -c gas_cwrap.d $(PIC) $(DFLAGS) -op \
		$(GAS_FILES) $(GAS_LUA_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) $(GASDYN_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(UTIL_FILES) $(EQC_SRC_FILES) $(GEOM_FILES) $(GZIP_FILES)
	$(DMD) $(OF)libgas.so gas_cwrap.o $(PIC) $(DFLAGS) -op -shared \
		$(GAS_FILES) $(GAS_LUA_FILES) $(KINETICS_FILES) $(GASDYN_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(UTIL_FILES) $(EQC_SRC_FILES) $(GEOM_FILES) $(GZIP_FILES) $(KINETICS_LUA_FILES) \
		$(LIBEQC) $(LIBLUA) \
		-L-L=$(INSTALL_DIR)/lib $(DLINKFLAGS)

# Priority:
#   1. VIRTUAL_ENV - Standard set by python virtual-environments
#   2. CONDA_PREFIX - Set by (ana)conda environments
#   3. INSTALL_DIR/.venv - Create our own if none already exist
CONDA_PREFIX ?= $(INSTALL_DIR)/.venv
VIRTUAL_ENV ?= $(CONDA_PREFIX)

# Override if using uv with `PIP='uv pip'`
PIP ?= $(VIRTUAL_ENV)/bin/pip
PYTHON ?= python3

# If it doesn't already exist
$(VIRTUAL_ENV):
	- mkdir -p $(INSTALL_DIR)
	$(info Creating virtual environment ...)
	$(PYTHON) -m venv $(VIRTUAL_ENV)
	$(info Updating pip)
	$(VIRTUAL_ENV)/bin/pip install --upgrade pip -q

$(PY_PROGRAMS): $(VIRTUAL_ENV)
	make -C $(GDTK_PY_DIR) build
	$(info Installing python tools)
	$(PIP) install .
	- mkdir -p $(INSTALL_DIR)/bin
	- ln -sf $(VIRTUAL_ENV)/bin/$@ $(INSTALL_DIR)/bin/$@

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

$(LIBEQC):
	cd $(EQC_DIR); make

species-database/species-database.lua :
	cd species-database; make species-database.lua

# Some demo targets. These are left over from the early days of development.

gas_model_demo: gas_model_demo.d $(GAS_FILES) $(EQC_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) -g $(DFLAGS) gas_model_demo.d $(OF)gas_model_demo \
		$(GAS_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

ideal_gas_demo: ideal_gas_demo.d $(GAS_FILES) $(EQC_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) -g $(DFLAGS) ideal_gas_demo.d $(OF)ideal_gas_demo \
		$(GAS_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

luagas_model_demo: luagas_model_demo.d $(GAS_FILES) $(EQC_SRC_FILES) $(GAS_LUA_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(LIBEQC) $(LIBLUA)
	$(DMD) -g $(DFLAGS) luagas_model_demo.d $(OF)luagas_model_demo \
		$(GAS_FILES) $(EQC_SRC_FILES) $(GAS_LUA_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) $(LIBEQC) $(LIBLUA) \
		$(DLINKFLAGS)

# Jonathan's programs that were moved out of nm.

build_tree: build_tree.d $(GAS_FILES) $(EQC_SRC_FILES) $(GAS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) build_tree.d -dip1008 $(OF)build_tree \
		$(GAS_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

univariate_lut_demo: univariate_lut_demo.d $(GAS_FILES) $(EQC_SRC_FILES) $(GAS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) -debug -dip1008 univariate_lut_demo.d $(OF)univariate_lut_demo \
		$(GAS_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

# Individual test targets

cea_therm_cond_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

chemkin_therm_cond_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

chemkin_viscosity_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

cea_viscosity_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

cea_thermo_curves_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

cea_thermo_curves_complex_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)cea_thermo_curves_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)


init_gas_model_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

init_gas_model_complex_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)init_gas_model_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

ideal_gas_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

ideal_gas_complex_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)ideal_gas_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

ideal_helium_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

cubic_gas_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

cea_gas_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

perf_gas_mix_eos_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

sutherland_therm_cond_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

sutherland_viscosity_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

therm_perf_gas_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

therm_perf_gas_complex_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)therm_perf_gas_test $(DVERSION)complex_numbers \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

therm_perf_gas_equil_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

therm_perf_gas_mix_eos_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

two_temperature_gas_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

multi_temperature_gas_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

two_temperature_trans_props_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

multi_temperature_trans_props_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

composite_gas_test : $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)


very_viscous_air_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

wilke_mixing_therm_cond_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

wilke_mixing_viscosity_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

co2gas_sw_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

uniform_lut_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

uniform_lut_plus_ideal_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

adaptive_lut_CEA_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

lut_comp_demo: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA) lut_comp_demo.d
	$(DMD) $(OF)$@ $(DFLAGS) lut_comp_demo.d \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

ideal_gas_ab_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

vib_specific_nitrogen_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

vib_specific_co_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

vib_specific_co_mixture_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

two_temperature_reacting_argon_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

two_temperature_argon_plus_ideal_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

ideal_dissociating_gas_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

fuel_air_mix_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

equilibrium_gas_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

steam_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

two_temperature_gasgiant_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

gasgiant_transport_properties_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

therm_perf_gas_mix_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

gas_mixtures_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ diffusion/transport_properties_model.d diffusion/gas_mixtures.d \
		$(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)


two_temp_air_ci_dump: $(GAS_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)



