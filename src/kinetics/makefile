# makefile for kinetics module
# Can be used to perform unit tests
# and build stand-alone programs.

DMD ?= ldc2

INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR ?= ../../build

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a

include ../common.mk
include $(wildcard *.deps)

build-prep-chem: prep_chem.lua reaction.lua lex_elems.lua
	$(INSTALL_DATA) reaction.lua lex_elems.lua -t $(BUILD_DIR)/lib/ 
	$(INSTALL_EXEC) prep_chem.lua $(BUILD_DIR)/bin/prep-chem

build-chemkin2eilmer: chemkin2eilmer.lua lex_elems.lua reaction.lua
	$(INSTALL_DATA) reaction.lua lex_elems.lua -t $(BUILD_DIR)/lib/
	$(INSTALL_EXEC) chemkin2eilmer.lua $(BUILD_DIR)/bin/chemkin2eilmer

build-prep-kinetics: prep_kinetics.lua mechanism.lua lex_elems.lua
	$(INSTALL_DATA) mechanism.lua lex_elems.lua -t $(BUILD_DIR)/lib/
	$(INSTALL_EXEC) prep_kinetics.lua $(BUILD_DIR)/bin/prep-kinetics

d_previews := in dip1008
d_versions :=

imports := eqc
import_paths := ../extern/eqc/source

imports += gas kinetics nm ntypes util
import_paths += ..

extra_deps := $(LIBLUA) $(LIBEQC)

test: test_runner
	@echo "Running real tests"
	- exec ./test_runner --module kinetics

test_runner: FLAVOUR = debug
test_runner: compiler_flags += -unittest
test_runner: ../util/test_runner.d package.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

clean:
	- rm -f *.o *.mod *.obj *.deps
	- rm -f test_runner
	- rm -f two_temperature_argon_kinetics_test_results.data
	- rm -rf $(BUILD_DIR)/*
	- $(MAKE) -C $(LUA_DIR) clean
	- $(MAKE) -C $(EQC_DIR) clean

$(LIBLUA):
	$(MAKE) -C $(LUA_DIR) guess local

$(LIBEQC):
	$(MAKE) -C $(EQC_DIR)

