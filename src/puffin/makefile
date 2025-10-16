# makefile for Puffin, a 2D supersonic flow calculator.
# PJ 2022-01-21

DMD ?= ldc2

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug
WITH_FPE ?= 1

PROGRAMS := puffin puffin-prep puffin-post lrkt-run lrkt-prep lrkt-post

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a

# The install destination.
INSTALL_DIR ?= $(HOME)/gdtkinst

# The build destination sits locally for present
BUILD_DIR := ../../build

# The build destination sits locally for present
REVISION_STRING := $(shell git rev-parse --short HEAD 2>/dev/null || echo "unknown")

default: $(PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Puffin code built."

install:
	make -C ../ install-puffin

clean:
	- rm -f *.o *.deps
	- rm -f $(PROGRAMS)
	- $(MAKE) -C $(LUA_DIR) clean;
	- $(MAKE) -C $(EQC_DIR) clean

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

string_imports := .

# Non-d files to include
extra_deps := $(LIBLUA) $(LIBEQC)

puffin: puffin_main.d $(extra_deps) | buildinfo.json
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

puffin-prep puffin-post: puffin-%: puffin_%.py
	$(INSTALL_EXEC) $< $@

lrkt-run: lorikeet_main.d $(extra_deps) | buildinfo.json
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

lrkt-prep lrkt-post: lrkt-%: lorikeet_%.py
	$(INSTALL_EXEC) $< $@

$(LIBLUA):
	$(MAKE) -C $(LUA_DIR) guess local

$(LIBEQC):
	$(MAKE) -C $(EQC_DIR)
