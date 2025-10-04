# makefile for Puffin, a 2D supersonic flow calculator.
# PJ 2022-01-21

DMD ?= ldc2

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug
WITH_FPE ?= 1

PROGRAMS := puffin puffin-prep puffin-post lrkt-run lrkt-prep lrkt-post

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

REPO_DIR := $(shell cd ../../; pwd)
ifneq ("$(wildcard $(REPO_DIR)/.git)","")
    # Repository exists
    REVISION_STRING := $(shell git rev-parse --short HEAD)
else
    REVISION_STRING := unknown
endif

default: $(PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Puffin code built."

install: $(PROGRAMS) prep-gas prep-chem chemkin2eilmer libgas
	- mkdir -p $(INSTALL_DIR)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share
	cp $(PROGRAMS) $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/bin/* $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/lib/lpeg.so $(BUILD_DIR)/lib
	cp -r ../lib/* $(BUILD_DIR)/lib
	cp $(NML_LUA_MODULES) $(BUILD_DIR)/lib
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

clean:
	- rm *.o *.deps
	- rm $(PROGRAMS)
	- rm -r $(BUILD_DIR)/*
	- rm *_main_with_rev_string.d
	- $(MAKE) -C $(GAS_DIR) clean
	- $(MAKE) -C $(KINETICS_DIR) clean
	- $(MAKE) -C $(LUA_DIR) clean;
	- $(MAKE) -C $(EQC_DIR) clean

$(LIBLUA):
	$(MAKE) -C $(LUA_DIR) guess local

$(LIBEQC):
	$(MAKE) -C $(EQC_DIR)

include ../common.mk
include $(wildcard *.deps)

d_previews := in dip1008
d_versions := flavour_$(FLAVOUR) \
	$(if $(filter 1,$(WITH_FPE)), enable_fp_exceptions) \

# External dependencies
imports := eqc gzip
import_paths := \
	../extern/eqc/source \
	../extern/gzip

# Internal dependencies
imports += puffin gas gasdyn geom grid kinetics nm ntypes util
import_paths += ..

ifeq ($(WITH_DVODE),1)
	d_versions += with_dvode
	linker_libs += gfortran
endif

# Non-d files to include
extra_deps := $(LIBLUA) $(LIBEQC)

%_with_rev_string.d: %.d
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		$< > $@

puffin: puffin_main_with_rev_string.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

puffin-prep puffin-post: puffin-%: puffin_%.py
	$(INSTALL_EXEC) $< $@

lrkt-run: lorikeet_main_with_rev_string.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

lrkt-prep lrkt-post: lrkt-%: lorikeet_%.py
	$(INSTALL_EXEC) $< $@

prep-gas libgas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) build-$@

chemkin2eilmer prep-chem:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-$@
