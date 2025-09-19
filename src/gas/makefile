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
build-libgas: libgas.so gas-calc build-uniform-lut.py $(LIBEQC)
	- mkdir -p $(BUILD_DIR)/bin
	cp gas-calc $(BUILD_DIR)/bin/
	cp build-uniform-lut.py $(BUILD_DIR)/bin/build-uniform-lut; chmod +x $(BUILD_DIR)/bin/build-uniform-lut
	cp ../../examples/estcn/estcn.py $(BUILD_DIR)/bin/estcn; chmod +x $(BUILD_DIR)/bin/estcn
	- mkdir -p $(BUILD_DIR)/lib
	cp libgas.so $(BUILD_DIR)/lib/
	cp $(EQC_DIR)/libeqc.so $(BUILD_DIR)/lib/
	cp -r ../lib/* $(BUILD_DIR)/lib/

install: build-prep-gas build-libgas
	- mkdir -p $(INSTALL_DIR)
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

.PHONY: test-all test-real
test-all: test-real

test-real: test_runner_real
	@echo "Running tests"
	- ./test_runner_real --module gas

# We need to pull the test_runner.d to the front to ensure no linker errors
test_runner_real: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
		$(KINETICS_FILES) $(KINETICS_LUA_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(UTIL_FILES) $(LIBEQC) $(LIBLUA)
	@echo "Compiling test runner..."
	$(DMD) $(OF)$@ $(DFLAGS) --unittest $(DLINKFLAGS) \
		$(filter %/test_runner.d,$^) $(filter-out %/test_runner.d,$^)

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm -f *.o *.so *.a *.obj *.lib
	- rm -f *.mod
	- rm -rf $(BUILD_DIR)/*
	- rm -rf __pycache__
	- rm -f ${DEMO_PROGRAMS}
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

univariate_lut_demo: univariate_lut_demo.d $(GAS_FILES) $(EQC_SRC_FILES) $(GAS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) -debug -dip1008 univariate_lut_demo.d $(OF)univariate_lut_demo \
		$(GAS_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

# Individual test targets

two_temp_air_ci_dump: $(GAS_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
		$(GAS_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

