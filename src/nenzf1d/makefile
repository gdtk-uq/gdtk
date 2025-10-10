# makefile for nenzf1d (mostly copied from the l1d makefile).
# PJ 2020-09-26

DMD ?= ldc2

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug
WITH_FPE ?= 1

PROGRAMS := nenzf1d

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

install: 
	make -C ../ install-nenzf1d

clean:
	- rm -f *.o *.so *.deps
	- rm -f $(PROGRAMS)
	- rm -f *_with_rev_string.d
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

$(LIBLUA):
	$(MAKE) -C $(LUA_DIR) guess local

$(LIBEQC):
	$(MAKE) -C $(EQC_DIR)
