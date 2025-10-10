# makefile for slf 
# NNG 2024-05-15

DMD ?= ldc2

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug
WITH_FPE ?= 0

PROGRAMS := slf libslf.so
SLF_FILES := io.d linalg.d misc.d

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a

GAS_DIR := ../gas
KINETICS_DIR := ../kinetics

DYAML_DIR := $(DGD_REPO)/src/extern/D-YAML/source/dyaml
TINYENDIAN_DIR := $(DGD_REPO)/src/extern/tinyendian/source

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

default: $(PROGRAMS) $(SLF_FILES)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "slf code built."

install: $(PROGRAMS)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/lib/gdtk
	- cp slf $(BUILD_DIR)/bin
	- cp libslf.so $(BUILD_DIR)/lib
	- cp slf.py $(BUILD_DIR)/lib/gdtk
	cp $(LUA_DIR)/install/bin/* $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/lib/lpeg.so $(BUILD_DIR)/lib
	cp -r ../lib/* $(BUILD_DIR)/lib
	@echo "Installing to $(INSTALL_DIR)"
	- mkdir -p $(INSTALL_DIR)
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

clean:
	- rm *.o *.so *deps
	- rm $(PROGRAMS)
	- rm -r $(BUILD_DIR)/*
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(EQC_DIR); make clean
	- cd $(DYAML_DIR); rm -rf *.o
	- cd $(TINYENDIAN_DIR); rm -rf *.o

include ../common.mk
include $(wildcard *.deps)

d_previews := in dip1008
d_versions := flavour_$(FLAVOUR) complex_numbers

ifeq ($(WITH_FPE),1)
    d_versions += enable_fp_exceptions
endif

imports := dyaml eqc gzip tinyendian
import_paths := \
	../extern/D-YAML/source \
	../extern/eqc/source \
	../extern/gzip \
	../extern/tinyendian/source 

imports += slf gas kinetics nm ntypes util
import_paths += ..

extra_deps := $(LIBLUA) $(LIBEQC)

slf: main.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

%.o: target_type = staticLibrary
%.o: %.d
	$(DMD) -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

libslf.so: linker_paths += $(INSTALL_DIR)/lib
libslf.so: target_type := dynamicLibrary
libslf.so: cwrap.o 
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

$(LIBLUA):
	cd $(LUA_DIR); make guess local

$(LIBEQC):
	cd $(EQC_DIR); make
