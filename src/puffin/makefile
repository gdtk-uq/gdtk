# makefile for Puffin, a 2D supersonic flow calculator.
# PJ 2022-01-21

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

PROGRAMS := puffin lrkt-run
PY_PROGRAMS := puffin-prep puffin-post lrkt-prep lrkt-post

PUFFIN_FILES := config.d marching_calc.d streamtube.d \
	cell.d face.d flow.d flux.d interp.d

LORIKEET_FILES := config.d transient_calc.d fluidblock.d \
	cell.d face.d flow.d flux.d interp.d

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
        DFLAGS := -profile -w -O2 -release -inline -boundscheck=off -version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O2 -release -inline -boundscheck=off -version=flavour_fast
    endif
    OF := -of
    DVERSION := -version=
    DLINKFLAGS :=
    DLINKFLAGS := $(DLINKFLAGS) -L-ldl
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
        DFLAGS := -fprofile-generate -g -w -O2 -release -enable-inlining -boundscheck=off -d-version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O2 -release -enable-inlining -boundscheck=off -d-version=flavour_fast
        ifeq ($(PLATFORM), linux)
	    # Full link-time optimization does not play nicely on macOS
	    DFLAGS += -flto=full
        endif
    endif
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

ifeq ($(WITH_DVODE),1)
    DFLAGS += $(DVERSION)with_dvode
    DLINKFLAGS += -L-lgfortran
endif

default: $(PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Puffin code built."

install: $(PROGRAMS) $(PY_PROGRAMS) prep-gas prep-chem chemkin2eilmer libgas
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
	- rm *.o
	- rm $(PROGRAMS)
	- rm -rf __pycache__/ *.egg-info/
	- rm -r $(BUILD_DIR)/*
	- rm puffin_main_with_rev_string.d lorikeet_main_with_rev_string.d
	- cd $(LUA_DIR); make clean;
	- cd $(GAS_DIR); make clean; rm libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(EQC_DIR); make clean

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

$(LIBEQC):
	cd $(EQC_DIR); make

puffin: puffin_main.d $(PUFFIN_FILES) \
	$(GEOM_FILES) \
	$(GAS_FILES) $(EQC_SRC_FILES) $(LIBEQC) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		puffin_main.d > puffin_main_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)puffin \
		puffin_main_with_rev_string.d \
		$(PUFFIN_FILES) \
		$(GEOM_FILES) \
		$(GAS_FILES) $(EQC_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBEQC) $(LIBLUA) \
		$(DLINKFLAGS)

lrkt-run: lorikeet_main.d $(LORIKEET_FILES) \
	$(GEOM_FILES) \
	$(GAS_FILES) $(EQC_SRC_FILES) $(LIBEQC) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		lorikeet_main.d > lorikeet_main_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)lrkt-run \
		lorikeet_main_with_rev_string.d \
		$(LORIKEET_FILES) \
		$(GEOM_FILES) \
		$(GAS_FILES) $(EQC_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBEQC) $(LIBLUA) \
		$(DLINKFLAGS)

# Priority:
#   1. VIRTUAL_ENV - Standard set by python virtual-environments
#   2. CONDA_PREFIX - Set by (ana)conda environments
#   3. INSTALL_DIR/.venv - Create our own if none already exist
CONDA_PREFIX ?= $(INSTALL_DIR)/.venv
VIRTUAL_ENV ?= $(CONDA_PREFIX)

# Override if using uv with `PIP='uv pip'`
PIP ?= $(VIRTUAL_ENV)/bin/pip
PYTHON ?= python3

# If it doesn't already exist
$(VIRTUAL_ENV):
	- mkdir -p $(INSTALL_DIR)
	$(info Creating virtual environment ...)
	$(PYTHON) -m venv $(VIRTUAL_ENV)
	$(info Updating pip)
	$(VIRTUAL_ENV)/bin/pip install --upgrade pip -q

$(PY_PROGRAMS): $(VIRTUAL_ENV)
	$(info Installing python tools)
	$(PIP) install -e . --find-links $(INSTALL_DIR)/share/python-packages/
	- mkdir -p $(INSTALL_DIR)/bin
	- ln -sf $(VIRTUAL_ENV)/bin/$@ $(INSTALL_DIR)/bin/$@

prep-gas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) build-prep-gas

libgas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) build-libgas

prep-chem:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-prep-chem

chemkin2eilmer:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-chemkin2eilmer
