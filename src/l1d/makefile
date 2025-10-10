# makefile for Lagrangian 1D flow solver.
# PJ 2020-04-05

DMD ?= ldc2

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug
WITH_FPE ?= 1

PROGRAMS := l1d4 l1d4-prep l1d4-xt-viewer

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
	@echo "L1d4 simulation code built."

install:
	make -C ../ install-l1d

clean:
	- rm -f *.o *.deps
	- rm -f $(PROGRAMS)
	- rm -f main_with_rev_string.d
	- $(MAKE) -C $(LUA_DIR) clean;
	- $(MAKE) -C $(EQC_DIR) clean;

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
imports += l1d gas gasdyn geom kinetics nm ntypes util
import_paths += ..

# Non-d files to include
extra_deps := $(LIBLUA) $(LIBEQC)

main_with_rev_string.d: main.d
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		main.d > main_with_rev_string.d

l1d4: main_with_rev_string.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

l1d4-prep: l1d4-prep.py
	$(INSTALL_EXEC) l1d4-prep.py l1d4-prep

l1d4-xt-viewer: xt-viewer.py
	$(INSTALL_EXEC) xt-viewer.py l1d4-xt-viewer

prep-gas libgas:
	$(MAKE) -C $(GAS_DIR) BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) build-$@

prep-chem prep-kinetics chemkin2eilmer:
	$(MAKE) -C $(KINETICS_DIR) BUILD_DIR=$(BUILD_DIR) build-$@

$(LIBLUA):
	$(MAKE) -C $(LUA_DIR) guess local

$(LIBEQC):
	$(MAKE) -C $(EQC_DIR)
