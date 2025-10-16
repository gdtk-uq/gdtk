# makefile for the gas module

DMD ?= ldc2

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug
WITH_FPE ?= 1

INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR := ../../build

DEMO_PROGRAMS := ideal_gas_demo luagas_model_demo
# FIXME: Neither of these demos compile due to changed source code.
# The files they rely on have changed or been removed.
# - gas_model_demo
# - univariate_lut_demo

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a

# ----------------------------------------------------------------------
# Here begins the list of targets, starting with the top-level actions.
# ----------------------------------------------------------------------

CIF_DIR := species-database/collision-integrals
COLLISION_INTEGRAL_FILES := $(CIF_DIR)/gupta_etal_1990_CI_data.lua \
	$(CIF_DIR)/wright_etal_CI_data.lua \
	$(CIF_DIR)/palmer_etal_CI_data.lua

DATA_FILES := $(COLLISION_INTEGRAL_FILES) \
	species-database/species-list.txt \
	species-database/species-database.lua

install: 
	make -C ../ install-gas

demo: $(DEMO_PROGRAMS)
	echo "Demo programs built."

# prep-gas and its database is required by the Eilmer flow solver.
build-prep-gas: prep_gas.lua species_data_converter.lua $(DATA_FILES)
	- mkdir -p $(BUILD_DIR)/bin $(BUILD_DIR)/data
	$(INSTALL_EXEC) prep_gas.lua $(BUILD_DIR)/bin/prep-gas
	$(INSTALL_EXEC) species_data_converter.lua $(BUILD_DIR)/bin/species-data-converter
	cp -r $(DATA_FILES) $(BUILD_DIR)/data/

# The following are not required by the Eilmer flow solver but
# they are useful for building stand-alone flow-analysis programs.
build-libgas: libgas.so gas-calc build-uniform-lut.py
	- mkdir -p $(BUILD_DIR)/bin $(BUILD_DIR)/lib
	$(INSTALL_EXEC) gas-calc $(BUILD_DIR)/bin/
	$(INSTALL_EXEC) build-uniform-lut.py $(BUILD_DIR)/bin/build-uniform-lut
	$(INSTALL_EXEC) ../../examples/estcn/estcn.py $(BUILD_DIR)/bin/estcn
	$(INSTALL_DATA) libgas.so $(BUILD_DIR)/lib/

clean:
	- rm -f *.o *.so *.a *.obj *.lib *.deps
	- rm -f *.mod
	- rm -rf __pycache__
	- rm -f $(DEMO_PROGRAMS)
	- rm -f gas-calc
	- rm -f test_runner_real
	- rm -f species-database/species-database.lua
	# Clean up test files
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
	- $(MAKE) -C $(LUA_DIR) clean
	- $(MAKE) -C $(EQC_DIR) clean

# Specific targets for individual artefacts.

include ../common.mk

d_previews := in dip1008
d_versions := flavour_$(FLAVOUR) \
	$(if $(filter 1,$(WITH_FPE)), enable_fp_exceptions)

# External dependencies
imports := eqc gzip
import_paths := \
	../extern/eqc/source \
	../extern/gzip

# Internal dependencies (all from gdtk/src)
imports += gas gasdyn geom grid kinetics nm ntypes util
import_paths += ..

# Non-d files to include
extra_deps := $(LIBLUA) $(LIBEQC)

# Include all possible target requirement files
include $(wildcard *.deps)

.PHONY: test-all test-real
test-all: test-real

test-real: gas-calc test_runner_real
	@echo "Running tests"
	- ./test_runner_real --module gas

# We need to pull the test_runner.d to the front to ensure no linker errors
test_runner_real: FLAVOUR = debug
test_runner_real: compiler_flags += -unittest
test_runner_real: test_gas_calc.d package.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

# compiler_flags is lazily defined in common.mk
gas-calc: gas_calc.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

%.o: target_type = staticLibrary
%.o: %.d
	$(DMD) -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

libgas.so: linker_paths += $(INSTALL_DIR)/lib
libgas.so: target_type = dynamicLibrary
libgas.so: gas_cwrap.o $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

species-database/species-database.lua :
	$(MAKE) -C ./species-database species-database.lua

prep-gas: prep_gas.lua
	$(INSTALL_EXEC) $< $@

species-data-converter: species_data_converter.lua
	$(INSTALL_EXEC) $< $@

build-uniform-lut: build-uniform-lut.py
	$(INSTALL_EXEC) $< $@

estcn: ../../examples/estcn/estcn.py
	$(INSTALL_EXEC) $< $@

# Some demo targets. These are left over from the early days of development.
# Since the target and filename are the same, these recipes are super simple!
$(DEMO_PROGRAMS): %: %.d $(extra_deps)
	$(DMD) -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

# Individual test targets
two_temp_air_ci_dump: d_versions += two_temp_air_ci_dump
two_temp_air_ci_dump: two_temperature_air.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

$(LIBLUA):
	$(MAKE) -C $(LUA_DIR) guess local

$(LIBEQC):
	$(MAKE) -C $(EQC_DIR)

