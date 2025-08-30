# makefile for Eilmer
# Builds main simulation program by default.
# For example:
# make install
#
# We need some trickery in make to get a newline into a string. We define it as a variable.
define nl


endef

DMD_ALLOWED_OPTS := dmd ldc2
# We can specify the DMD64 compiler as DMD=dmd on the command-line
# when invoking this makefile.
DMD ?= ldc2
ifeq ($(filter $(DMD), $(DMD_ALLOWED_OPTS)),)
    $(warning DMD=$(DMD) is not available as an option for setting the D compiler.)
    $(warning Available options are:${nl}    $(DMD_ALLOWED_OPTS))
    $(error Exiting unsuccessfully.${nl}Please fix DMD selection on command line or as environment variable setting.)
endif
NVCC ?= nvcc
GPP ?= g++
# default is NOT to opt out of MPI build
OPT_OUT_MPI ?= 0
MPI_IMPLEMENTATION := OpenMPI
# default is NOT to use MPICH
WITH_MPICH ?= 0

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
# As default, we compile with the 'debug' option which produces a code
# that will have better checks and more detailed error messages.
# The 'fast' option requests production mode and to include optimisations.
FLAVOUR_ALLOWED_OPTS := debug fast profile
FLAVOUR ?= debug
ifeq ($(filter $(FLAVOUR), $(FLAVOUR_ALLOWED_OPTS)),)
    $(warning FLAVOUR=$(FLAVOUR) is not available as an option for building.)
    $(warning Available options are: ${nl}    $(FLAVOUR_ALLOWED_OPTS))
    $(error Exiting unsuccessfully.${nl}Please fix FLAVOUR selection.)
endif
WITH_FPE ?= 1

# PLATFORM options are linux, macosx
ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

WITH_OPENCL_GPU_CHEM ?= 0
WITH_CUDA_GPU_CHEM ?= 0
DEBUG_CHEM ?= 0
WITH_DIAGNOSTICS ?= 0
MULTI_SPECIES_GAS ?= 1
MULTI_T_GAS ?= 1
MHD ?= 1
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
AUX_PROGRAMS_SHARE += lmr-complete.sh gdtk-module

MPI_LIBRARY_DIRS = $(shell mpicc --showme:libdirs)
MPI_LIB_DIRS_SEARCH = $(foreach d, $(MPI_LIBRARY_DIRS), -L-L$d)

include lmr-files.mk

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

NML_DIR := ../lib
include $(NML_DIR)/nml_files.mk

GAS_DIR := ../gas
GAS_PACKAGE := gas
include $(GAS_DIR)/gas_files.mk

GRID_DIR := ../grid_utils
include $(GRID_DIR)/grid_utils_files.mk

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

GEOM_DIR := ../geom
include $(GEOM_DIR)/geom_files.mk

GASDYN_DIR := ../gasdyn
include $(GASDYN_DIR)/gasdyn_files.mk

EQC_DIR := ../extern/eqc/source
include $(EQC_DIR)/eqc_files.mk

NTYPES_DIR := ../ntypes
NTYPES_FILES := $(NTYPES_DIR)/complex.d

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

DYAML_DIR := ../extern/D-YAML/source/dyaml
include $(DYAML_DIR)/dyaml_files.mk

TINYENDIAN_DIR := ../extern/tinyendian/source
include $(TINYENDIAN_DIR)/tinyendian_files.mk

MPL_DIR := ../extern/matplotlib.d/source
MPL_FILES := $(MPL_DIR)/matplotlibd/pyplot.d \
	$(MPL_DIR)/matplotlibd/core/pycall.d \
	$(MPL_DIR)/matplotlibd/core/translate.d

GPERF_DIR := ../extern/gperftools_d/source/gperftools_d
GPERF_FILES := $(GPERF_DIR)/heap_profiler.d \
	$(GPERF_DIR)/malloc_extension_c.d \
	$(GPERF_DIR)/malloc_hook_c.d \
	$(GPERF_DIR)/profiler.d \
	$(GPERF_DIR)/stacktrace.d \
	$(GPERF_DIR)/tcmalloc.d \

# The install destination.
INSTALL_DIR ?= $(HOME)/gdtkinst

# The build destination sits locally for present
BUILD_DIR := ./build
BUILD_DATE := $(shell date)

REVISION_STRING := $(shell git rev-parse --short HEAD)
FULL_REVISION_STRING := $(shell git rev-parse HEAD)
REVISION_AGE := $(shell git log -1 --format=%cd --date=relative)
REVISION_DATE := $(shell git log -1 --format=%cd)
REPO_DIR := $(shell cd ../../; pwd)

DFLAGS :=
DLINKFLAGS :=

ifeq ($(DMD), dmd)
    DEBUG_DFLAGS := -w -g -debug -version=flavour_debug
    PROFILE_DFLAGS := -profile -w -g -O -release -boundscheck=off -version=flavour_profile
    FAST_DFLAGS := -w -g -O -release -boundscheck=off -version=flavour_fast
    OF := -of
    DVERSION := -version=
    ifeq ($(findstring with_libplot,$(LIBPLOT_VERSION_STR)), with_libplot)
        DFLAGS := $(DFLAGS) $(DVERSION)$(LIBPLOT_VERSION_STR)
        DLINKFLAGS := $(DLINKFLAGS) -L-lplot
    endif
	DLINKFLAGS := $(DLINKFLAGS) -L-ldl
	ifeq ($(PLATFORM), linux)
		# For note on "-Wl,-E", see: https://web.tecgraf.puc-rio.br/~lhf/ftp/lua/install.html#adding
		# DLINKFLAGS += -L-Wl,-E
		# 2025-04-05 Linker does not recognise this flag on Linux Mint 22 when using the dmd compiler.
	endif
    ifeq ($(PLATFORM), macosx)
        DLINKFLAGS += -L-ld_classic
    endif
endif
ifeq ($(DMD), ldc2)
    # -fprofile-generate will result in profraw files being written
    # that may be viewed, showing the top 10 functions with internal block counts
    # llvm-profdata show -text -topn=10 <profraw-file>
    DEBUG_DFLAGS := -w -g --d-debug --d-version=flavour_debug
    PROFILE_DFLAGS := -fprofile-generate -g -w -O2 -release -enable-inlining -boundscheck=off --d-version=flavour_profile
    ifeq ($(PLATFORM), macosx)
        # 2025-01-15: -O2 messes with the piston-in-tube example on MacOS so we back off to -O1
        FAST_DFLAGS := -w -g -O1 -release -enable-inlining -boundscheck=off --d-version=flavour_fast
    else
        # The usual FAST flags for linux.
        FAST_DFLAGS := -w -g -O2 -release -enable-inlining -boundscheck=off --d-version=flavour_fast
    endif
    ifeq ($(PLATFORM), linux)
	  # Full link-time optimization does not play nicely on macOS
	  FAST_DFLAGS += -flto=thin
    endif
    ifeq ($(WITH_THREAD_SANITIZER), 1)
        DFLAGS := $(DFLAGS) -fsanitize=thread
    endif
    ifeq ($(WITH_ADDRESS_SANITIZER), 1)
        DFLAGS := $(DFLAGS) -fsanitize=address
    endif
    OF := -of=
    DVERSION := -d-version=
    DLINKFLAGS :=
    ifeq ($(WITH_THREAD_SANITIZER), 1)
        DLINKFLAGS := $(DLINKFLAGS) -fsanitize=thread
    endif
    ifeq ($(WITH_ADDRESS_SANITIZER), 1)
        DLINKFLAGS := $(DLINKFLAGS) -fsanitize=address
    endif
    #ifeq ($(FLAVOUR), profile)
    #    DLINKFLAGS := $(DLINKFLAGS) -Wl,-fprofile-generate
    #endif
    ifeq ($(findstring with_libplot,$(LIBPLOT_VERSION_STR)), with_libplot)
        DFLAGS := $(DFLAGS) $(DVERSION)$(LIBPLOT_VERSION_STR)
        DLINKFLAGS := $(DLINKFLAGS) -L-lplot
    endif
	DLINKFLAGS := $(DLINKFLAGS) -L-ldl
	ifeq ($(PLATFORM), linux)
		# For note on "-Wl,-E", see: https://web.tecgraf.puc-rio.br/~lhf/ftp/lua/install.html#adding
		DLINKFLAGS += -L-Wl,-E
	endif
    ifeq ($(PLATFORM), macosx)
		DLINKFLAGS += -L-ld_classic
    endif
endif

# Set the flavour to be compiler correct flags
ifeq ($(FLAVOUR), debug)
    FLAVOUR_FLAGS := $(DEBUG_DFLAGS)
endif
ifeq ($(FLAVOUR), profile)
    FLAVOUR_FLAGS := $(PROFILE_DFLAGS)
endif
ifeq ($(FLAVOUR), fast)
    FLAVOUR_FLAGS := $(FAST_DFLAGS)
    WITH_FPE ?= 0
endif

# Add platform version if necessary
ifeq ($(PLATFORM), macosx)
	DFLAGS += $(DVERSION)macosx
endif

# DIP1008 allows throwing of exceptions in @nogc code. Appeared in 2.079.0, 2018-03-01.
# This rules out the use of gdc for compiling the code.
DFLAGS += -dip1008

DFLAGS += -I.. \
		  -I$(GAS_DIR) \
		  -I$(GZIP_DIR) \
		  -I$(GRID_DIR) \
		  -I$(GEOM_DIR) \
		  -I$(UTIL_DIR) \
		  -I$(NM_DIR) \
		  -I$(NTYPES_DIR) \
		  -I$(KINETICS_DIR) \
		  -I$(GASDYN_DIR) \
		  -I$(EQC_DIR) \
		  -I$(DYAML_DIR)/.. \
		  -I$(TINYENDIAN_DIR)

ifeq ($(DEBUG_CHEM),1)
    DFLAGS += $(DVERSION)debug_chem
endif

ifeq ($(WITH_FPE),1)
    DFLAGS += $(DVERSION)enable_fp_exceptions
endif

ifeq ($(MULTI_SPECIES_GAS),1)
    DFLAGS += $(DVERSION)multi_species_gas
endif

ifeq ($(MULTI_T_GAS),1)
    DFLAGS += $(DVERSION)multi_T_gas
endif

ifeq ($(MHD),1)
    DFLAGS += $(DVERSION)MHD
endif

ifeq ($(TURBULENCE),1)
    DFLAGS += $(DVERSION)turbulence
endif

DLINKFLAGS += $(LINKER_FLAGS)

include $(LMR_DIR)lmr-rules.mk
