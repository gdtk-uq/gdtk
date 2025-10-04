# makefile for Eilmer
# Builds main simulation program by default.
# For example:
# make install
#
# We need some trickery in make to get a newline into a string. We define it as a variable.
define nl


endef

NVCC ?= nvcc
GPP ?= g++
# default is NOT to opt out of MPI build
OPT_OUT_MPI ?= 0
MPI_IMPLEMENTATION := OpenMPI
# default is NOT to use MPICH
WITH_MPICH ?= 0

WITH_FPE ?= 1

# PLATFORM options are linux, macosx
ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
    USING_MAC := 1
else
    PLATFORM := linux
    USING_MAC :=
endif
$(info PLATFORM=$(PLATFORM))

WITH_OPENCL_GPU_CHEM ?= 0
WITH_CUDA_GPU_CHEM ?= 0
DEBUG_CHEM ?= 0
WITH_DIAGNOSTICS ?= 0
MULTI_SPECIES_GAS ?= 1
MULTI_T_GAS ?= 1
MHD ?= 1
NK ?= 1
TURBULENCE ?= 1
WITH_E4DEBUG ?= 0
WITH_CHECK_JAC ?=0

# For the ldc2 compiler, we can turn on some extra checking at run-time.
# There will be a computational cost associated with these checks.
WITH_THREAD_SANITIZER ?= 0
WITH_ADDRESS_SANITIZER ?= 0

OPENMPI_DIR := ../extern/OpenMPI
OPENMPI_FILES := $(OPENMPI_DIR)/source/mpi/package.d

MPI_DIR = $(OPENMPI_DIR)
MPI_FILES = $(OPENMPI_FILES)

MPICH_DIR := ../extern/cray-mpich
MPICH_FILES := $(MPICH_DIR)/mpi.d

GAS_DIR = ../gas
GEOM_DIR = ../geom
GRID_DIR = ../grid_utils
KINETICS_DIR = ../kinetics
NML_DIR := ../lib
include $(NML_DIR)/nml_files.mk

PROGRAMS := lmr lmr-debug
SUB_PROGRAMS := lmr-run lmrZ-run
ifeq ($(OPT_OUT_MPI), 0)
	SUB_PROGRAMS += lmr-mpi-run lmrZ-mpi-run
endif
ifeq ($(WITH_MPICH),1)
    MPI_LIBRARY_DIRS = $(shell mpicc -link_info | cut -d' ' -f 3)
    MPI_LIB_DIRS_SEARCH = $(foreach d, $(MPI_LIBRARY_DIRS), -L$d)
    MPI_DIR = $(MPICH_DIR)
    MPI_FILES = $(MPICH_FILES)
	MPI_IMPLEMENTATION := MPICH
endif
ifeq ($(WITH_CHECK_JAC), 1)
	SUB_PROGRAMS += lmr-check-jacobian lmrZ-check-jacobian
endif
SHARE_FILES :=
ETC_FILES := lmr.cfg
PY_PROG_DIR := python-programs
PY_PROGRAMS := lmr-verify
AUX_PROGRAMS := prep-gas ugrid_partition prep-chem chemkin2eilmer prep-kinetics
AUX_PROGRAMS += gdtk-module

MPI_LIBRARY_DIRS = $(shell mpicc --showme:libdirs)
MPI_LIB_DIRS_SEARCH = $(foreach d, $(MPI_LIBRARY_DIRS), -L-L$d)

include lmr-files.mk

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a
include $(EQC_DIR)/eqc_files.mk

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/install/lib

# The install destination.
INSTALL_DIR ?= $(HOME)/gdtkinst

# The build destination sits locally for present
BUILD_DIR := ../../build
BUILD_DATE := $(shell date)

REVISION_STRING := $(shell git rev-parse --short HEAD)
FULL_REVISION_STRING := $(shell git rev-parse HEAD)
REVISION_AGE := $(shell git log -1 --format=%cd --date=relative)
REVISION_DATE := $(shell git log -1 --format=%cd)
REPO_DIR := $(shell cd ../../; pwd)

default: $(PROGRAMS) $(SUB_PROGRAMS)

# The pipe symbol (|) separates normal prerequisites from order-only ones.
# Order-only prerequisites being updated do *not* cause the target to
# become out of date. This prevents unneccesary rebuilds of targets when
# only the build-time changes in `buildinfo.json`.
#
# Including `buildinfo.json` as a order-only prerequisite here ensures that
# all programs have the same build information (i.e. build-time) when
# calling `make install`.
install: $(PROGRAMS) $(SUB_PROGRAMS) $(PY_PROGRAMS) $(AUX_PROGRAMS) | buildinfo.json
	- mkdir -p $(INSTALL_DIR)/bin
	- mkdir -p $(INSTALL_DIR)/lib
	- mkdir -p $(INSTALL_DIR)/etc
	- mkdir -p $(INSTALL_DIR)/share
	cp $(PROGRAMS) $(INSTALL_DIR)/bin
	cp $(SUB_PROGRAMS) $(INSTALL_DIR)/bin
	cp $(PY_PROGRAMS) $(INSTALL_DIR)/bin
	cp lua-modules/*.lua $(INSTALL_DIR)/lib/
	cp lmr.cfg $(INSTALL_DIR)/etc/
	cp $(LUA_DIR)/install/bin/* $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/lib/lpeg.so $(BUILD_DIR)/lib
	cp -r ../lib/* $(BUILD_DIR)/lib
	cp $(NML_LUA_MODULES) $(BUILD_DIR)/lib
	cp gdtk-module $(INSTALL_DIR)/share
	cp share/* $(INSTALL_DIR)/share
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

clean:
	- rm *.o
	- rm *.deps
	- rm -rf $(BUILD_DIR)/*
	- rm $(PROGRAMS)
	- rm $(SUB_PROGRAMS)
	- rm $(PY_PROGRAMS)
	- rm lmrconfig_with_str_subst.d
	- rm runsim_shared.d runsim_mpi.d
	- rm -f gdtk-module
	- rm -f lmr-check-jacobian lmrZ-check-jacobian
	- $(MAKE) -C $(OPENMPI_DIR) clean
	- $(MAKE) -C $(MPICH_DIR) clean
	- $(MAKE) -C $(GEOM_DIR) clean
	- $(MAKE) -C $(GZIP_DIR) clean
	- $(MAKE) -C $(GAS_DIR) clean
	- $(MAKE) -C $(KINETICS_DIR) clean
	- $(MAKE) -C $(GRID_DIR) clean
	- $(MAKE) -C $(LUA_DIR) clean
	- $(MAKE) -C $(EQC_DIR) clean

$(LIBGASF):
	$(MAKE) -C $(GAS_DIR) BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) libgasf.a

$(LIBLUA):
	$(MAKE) -C $(LUA_DIR) guess local

$(LIBEQC):
	$(MAKE) -C $(EQC_DIR)

$(MPI_FILES):
	$(MAKE) -C $(MPI_DIR)

# Makefiles don't like multiline strings. Using define and nl is
# the cleanest way to make a multiline variable.
define buildinfo
{
    "buildDate": "$(BUILD_DATE)",
    "revisionDate": "$(REVISION_DATE)",
    "revisionId": "$(REVISION_STRING)",
    "fullRevisionId": "$(FULL_REVISION_STRING)",
    "parallelType": "$(MPI_IMPLEMENTATION)"
}
endef

buildinfo.json:
	@echo "Creating 'buildinfo.json'..."
	# Ensure the define becomes a proper multi-line string.
	@printf '%b' '$(subst ${nl},\n,${buildinfo})' > $@

include ../common.mk
include $(wildcard *.deps)

d_previews := in dip1008
d_versions := flavour_$(FLAVOUR) \
	$(if $(USING_MAC),macosx) \
	$(if $(filter 1,$(DEBUG_CHEM)),debug_chem) \
	$(if $(filter 1,$(WITH_FPE)),enable_fp_exceptions) \
	$(if $(filter 1,$(MULTI_SPECIES_GAS)),multi_species_gas) \
	$(if $(filter 1,$(MULTI_T_GAS)),multi_T_gas) \
	$(if $(filter 1,$(MHD)),MHD) \
	$(if $(filter 1,$(NK)),newton_krylov) \
	$(if $(filter 1,$(TURBULENCE)),turbulence)

string_imports := .

# External dependencies
imports := dyaml eqc gzip tinyendian
import_paths := \
	../extern/D-YAML/source \
	../extern/eqc/source \
	../extern/gzip \
	../extern/tinyendian/source

ifeq ($(findstring with_libplot,$(LIBPLOT_VERSION_STR)), with_libplot)
    d_versions += $(LIBPLOT_VERSION_STR)
    linker_libs += plot
endif

# Internal dependencies
imports += lmr gas gasdyn geom grid kinetics nm ntypes util
import_paths += ..

# Non-d files to include
extra_deps = $(LIBLUA) $(LIBEQC)

# ============================================================================
#                      BUILD RECIPES FOR LMR VARIANTS
# ============================================================================ 

# List all our targets explicitly
lmr lmr-debug lmr-run lmrZ-run lmr-mpi-run lmrZ-mpi-run:
.SECONDEXPANSION:

# compiler_flags is lazily defined in common.mk
lmr-debug: FLAVOUR = debug
lmr-debug lmr: main.d $(extra_deps) | buildinfo.json
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

# Use the pattern rule % to match *all* lmrZ[...] targets.
# Add complex specialisations
lmrZ%run: d_versions += complex_numbers

# Add MPI specialisations
lmr%mpi-run: d_versions   += mpi_parallel
lmr%mpi-run: imports      += mpi
lmr%mpi-run: import_paths += ../extern/OpenMPI/source
lmr%mpi-run: linker_libs  += mpi
lmr%mpi-run: linker_paths += $(MPI_LIBRARY_DIRS)
lmr%mpi-run: extra_deps   += $(MPI_FILES)

lmr%run: d_versions += run_main
lmr%run: commands/runsim.d $$(extra_deps) | buildinfo.json
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

lmr%check-jacobian: commands/checkjacobian.d $(extra_deps) | buildinfo.json
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

#### Python programs
lmr-verify: $(PY_PROG_DIR)/lmr_verify.py
	cp $< $@
	chmod +x $@

#### Auxiliary programs but part of making an eilmer install useful
prep-gas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) PLATFORM=$(PLATFORM) build-prep-gas

ugrid_partition:
	- mkdir -p $(BUILD_DIR)/bin
	cd $(GRID_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) PLATFORM=$(PLATFORM) ugrid_partition
	cd $(GRID_DIR); cp ugrid_partition $(BUILD_DIR)/bin

prep-chem chemkin2eilmer prep-kinetics:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) PLATFORM=$(PLATFORM) build-$@

gdtk-module:
	sed -e 's+PUT_REVISION_STRING_HERE+$(REVISION_STRING)+' \
		-e 's+PUT_COMPILER_NAME_HERE+$(DMD)+' \
		-e 's+PUT_INSTALL_DIR_HERE+$(INSTALL_DIR)+' \
		-e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
		-e 's+PUT_BUILD_DATE_HERE+$(BUILD_DATE)+' \
		-e 's+PUT_REVISION_AGE_HERE+$(REVISION_AGE)+' \
		-e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
		./gdtk-module-template > gdtk-module

# Special Targets
# ===============
# The .INTERMEDIATE target will ensure that generated files are removed after
# a successful recipe *if* the files did not exist prior to calling the recipe.
.INTERMEDIATE: buildinfo.json
# The .PHONY target will ensure that recipes are *always* run, regardless of
# if their dependencies (if any) are out of date.
.PHONY: install clean
