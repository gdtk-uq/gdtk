# makefile for nenzf1d (mostly copied from the l1d makefile).
# PJ 2020-09-26

DMD ?= ldc2

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug
WITH_FPE ?= 1

PROGRAMS := nenzf1d

NML_DIR := ../lib
include $(NML_DIR)/nml_files.mk

GAS_DIR := ../gas
KINETICS_DIR := ../kinetics

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a

# The install destination.
INSTALL_DIR ?= $(HOME)/gdtkinst

# The build destination sits locally for present
BUILD_DIR := ../../build
BUILD_DATE := $(shell date)
REPO_DIR := $(shell cd ../../; pwd)
ifneq ("$(wildcard $(REPO_DIR)/.git)","")
    # Repository exists
    REVISION_STRING := $(shell git rev-parse --short HEAD)
    REVISION_DATE := $(shell git log -1 --format=%cd)
else
    REVISION_STRING := unknown
    REVISION_DATE := unknown
endif

default: $(PROGRAMS) libnenzf1d.so
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "nenzf1d code built."

install: $(PROGRAMS) libnenzf1d.so prep-gas prep-chem chemkin2eilmer prep-kinetics
	- mkdir -p $(INSTALL_DIR)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share
	cp $(PROGRAMS) $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/bin/* $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/lib/lpeg.so $(BUILD_DIR)/lib
	cp -r ../lib/* $(BUILD_DIR)/lib
	cp $(NML_LUA_MODULES) $(BUILD_DIR)/lib
	cp nenzf1d.py $(BUILD_DIR)/lib
	cp libnenzf1d.so $(BUILD_DIR)/lib 2>/dev/null || :
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

clean:
	- rm *.o *.so *.deps
	- rm $(PROGRAMS)
	- rm -r $(BUILD_DIR)/*
	- rm *_with_rev_string.d
	- $(MAKE) -C $(GAS_DIR) clean
	- $(MAKE) -C $(KINETICS_DIR) clean
	- $(MAKE) -C $(LUA_DIR) clean
	- $(MAKE) -C $(EQC_DIR) clean

include ../common.mk
include $(wildcard *.deps)

d_previews := in dip1008
d_versions := flavour_$(FLAVOUR) \
	$(if $(filter 1,$(WITH_FPE)), enable_fp_exceptions) \

# External dependencies
imports := dyaml eqc gzip tinyendian
import_paths := \
	../extern/D-YAML/source \
	../extern/eqc/source \
	../extern/gzip \
	../extern/tinyendian/source

# Internal dependencies
imports += nenzf1d gas gasdyn geom kinetics nm ntypes util
import_paths += ..

# Non-d files to include
extra_deps := $(LIBLUA) $(LIBEQC)

%_with_rev_string.d: %.d
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		$< > $@

nenzf1d: main_with_rev_string.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^

%.o: target_type = staticLibrary
%.o: %.d
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

# .so files generally don't have any deps other than their
# .o files, so dep tracking could be toggled off?
libnenzf1d.so: target_type = dynamicLibrary
libnenzf1d.so: linker_paths += $(INSTALL_DIR)/lib
libnenzf1d.so: cwrap.o $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

prep-gas:
	$(MAKE) -C $(GAS_DIR) BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) build-$@

prep-chem prep-kinetics chemkin2eilmer:
	$(MAKE) -C $(KINETICS_DIR) BUILD_DIR=$(BUILD_DIR) build-$@

$(LIBLUA):
	$(MAKE) -C $(LUA_DIR) guess local

$(LIBEQC):
	$(MAKE) -C $(EQC_DIR)
