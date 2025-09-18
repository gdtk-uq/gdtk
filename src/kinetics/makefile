# makefile for kinetics module
# Can be used to perform unit tests
# and build stand-alone programs.

DMD := ldmd2

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR ?= ../../build

NTYPES_DIR := ../ntypes
NTYPES_FILES := $(NTYPES_DIR)/complex.d

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk
LIBGASF := $(GAS_DIR)/libgasf.a

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a
include $(EQC_DIR)/eqc_files.mk

include kinetics_files.mk

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := ${LUA_DIR}/install/lib/liblua.a
LIBLUAPATH := ${LUA_DIR}/lib

DLINKFLAGS :=  -w -L-ldl
ifeq ($(PLATFORM), macosx)
    DLINKFLAGS += -L-ld_classic
endif

build-prep-chem: prep_chem.lua reaction.lua lex_elems.lua
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	cp reaction.lua lex_elems.lua $(BUILD_DIR)/lib/
	cp prep_chem.lua $(BUILD_DIR)/bin/prep-chem; chmod +x $(BUILD_DIR)/bin/prep-chem

build-chemkin2eilmer: chemkin2eilmer.lua lex_elems.lua reaction.lua
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	cp reaction.lua lex_elems.lua $(BUILD_DIR)/lib/
	cp chemkin2eilmer.lua $(BUILD_DIR)/bin/chemkin2eilmer; chmod +x $(BUILD_DIR)/bin/chemkin2eilmer

build-prep-kinetics: prep_kinetics.lua mechanism.lua lex_elems.lua
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	cp mechanism.lua lex_elems.lua $(BUILD_DIR)/lib/
	cp prep_kinetics.lua $(BUILD_DIR)/bin/prep-kinetics; chmod +x $(BUILD_DIR)/bin/prep-kinetics

test: test_runner
	@echo "Running real tests"
	- exec ./test_runner --module kinetics

# We need to pull the test_runner.d to the front to ensure no linker errors
test_runner: $(KINETICS_FILES) $(KINETICS_LUA_FILES) $(GAS_LUA_FILES) $(NM_FILES) $(NTYPES_FILES) $(EQC_SRC_FILES) \
		$(LIBEQC) $(LIBLUA) $(GAS_FILES) $(UTIL_FILES)
	@echo "Compiling test runner..."
	$(DMD) -of$@ -debug -g -w -dip1008 -unittest $(DLINKFLAGS) \
		$(filter %/test_runner.d,$^) $(filter-out %/test_runner.d,$^)

clean:
	- rm -f *.o *.mod *.obj
	- rm -f test_runner
	- rm -f two_temperature_argon_kinetics_test_results.data
	- rm -rf $(BUILD_DIR)/*
	- cd $(LUA_DIR); make clean
	- cd $(EQC_DIR); make clean

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

$(LIBEQC):
	cd $(EQC_DIR); make

