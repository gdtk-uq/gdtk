# makefile for nenzf1d (mostly copied from the l1d makefile).
# PJ 2020-09-26

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

PROGRAMS := nenzf1d

NENZF1D_FILES := configuration.d shock_tube_nozzle.d

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

NML_DIR := ../lib
include $(NML_DIR)/nml_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

GEOM_DIR := ../geom
include $(GEOM_DIR)/geom_files.mk

GASDYN_DIR := ../gasdyn
include $(GASDYN_DIR)/gasdyn_files.mk

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a
include $(EQC_DIR)/eqc_files.mk

NTYPES_DIR := ../ntypes
NTYPES_FILES := $(NTYPES_DIR)/complex.d

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib

DYAML_DIR := $(DGD_REPO)/src/extern/D-YAML/source/dyaml
include $(DYAML_DIR)/dyaml_files.mk

TINYENDIAN_DIR := $(DGD_REPO)/src/extern/tinyendian/source
include $(TINYENDIAN_DIR)/tinyendian_files.mk

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
    OF := -of
    DVERSION := -version=
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
        DFLAGS := -profile -w -O -release -inline -boundscheck=off -version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O -release -inline -boundscheck=off -version=flavour_fast
    endif
    PIC := -fPIC
    OF := -of
    DVERSION := -version=
    DLINKFLAGS := -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif
ifeq ($(DMD), ldc2)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -w -g -d-debug -d-version=flavour_debug
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
    OF := -of=
    DVERSION := -d-version=
    DLINKFLAGS := -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif
# DIP1008 allows throwing of exceptions in @nogc code. Appeared in 2.079.0, 2018-03-01.
# This rules out the use of gdc for compiling the code.
DFLAGS += -dip1008 -preview=in
DFLAGS += -I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GRID_DIR) -I$(GZIP_DIR)

ifeq ($(WITH_FPE),1)
    DFLAGS += $(DVERSION)enable_fp_exceptions
endif

default: $(PROGRAMS) libnenzf1d.so $(NENZF1D_FILES)
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
	- rm *.o *.so
	- rm $(PROGRAMS)
	- rm -r $(BUILD_DIR)/*
	- rm main_with_rev_string.d
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(EQC_DIR); make clean
	- cd $(DYAML_DIR); rm -rf *.o
	- cd $(TINYENDIAN_DIR); rm -rf *.o

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

$(LIBGASF):
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) libgasf.a

$(LIBEQC):
	cd $(EQC_DIR); make

nenzf1d: main.d $(NENZF1D_FILES) \
	$(GEOM_FILES) $(DYAML_FILES) $(TINYENDIAN_FILES) \
	$(GAS_FILES) $(EQC_SRC_FILES) $(LIBEQC) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		main.d > main_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)nenzf1d \
		main_with_rev_string.d \
		$(NENZF1D_FILES) $(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(GEOM_FILES) \
		$(GAS_FILES) $(EQC_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBEQC) $(LIBLUA) \
		$(DLINKFLAGS)

libnenzf1d.so: cwrap.d $(NENZF1D_FILES) \
	$(GEOM_FILES) $(DYAML_FILES) $(TINYENDIAN_FILES) \
	$(GAS_FILES) $(EQC_SRC_FILES) $(LIBEQC) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	$(DMD) -c cwrap.d $(PIC) $(DFLAGS) -op \
		$(NENZF1D_FILES) $(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(GEOM_FILES) \
		$(GAS_FILES) $(EQC_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBEQC) $(LIBLUA)
	$(DMD) $(OF)libnenzf1d.so cwrap.o $(PIC) $(DFLAGS) -op -shared \
		$(NENZF1D_FILES) $(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(GEOM_FILES) \
		$(GAS_FILES) $(EQC_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBEQC) $(LIBLUA) \
		-L-L=$(INSTALL_DIR)/lib $(DLINKFLAGS)

prep-gas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) build-prep-gas

prep-chem:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-prep-chem

chemkin2eilmer:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-chemkin2eilmer

prep-kinetics:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-prep-kinetics
