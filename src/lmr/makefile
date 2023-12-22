# makefile for Eilmer
# Builds main simulation program by default.
# For example:
# make install
#
# We can specify the DMD64 compiler as DMD=dmd on the command-line
# when invoking this makefile.
DMD ?= ldc2
NVCC ?= nvcc
GPP ?= g++

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
# As default, we compile with the 'debug' option which produces a code
# that will have better checks and more detailed error messages.
# The 'fast' option requests production mode and to include optimisations.
FLAVOUR ?= debug
WITH_FPE ?= 1

# PLATFORM options are linux, macosx
ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

# PARALLEL_FLAVOUR options are: shared, openmpi, mpich
PARALLEL_FLAVOUR ?= shared

WITH_MPI ?= 0
WITH_NK ?= 0
WITH_PIR ?= 0
WITH_SSC ?= 0
WITH_OPENCL_GPU_CHEM ?= 0
WITH_CUDA_GPU_CHEM ?= 0
DEBUG_CHEM ?= 0
WITH_COMPLEX_NUMBERS ?= 0
WITH_FPE ?= 0
WITH_DVODE ?= 0
WITH_DIAGNOSTICS ?= 0
MULTI_SPECIES_GAS ?= 1
MULTI_T_GAS ?= 1
MHD ?= 1
TURBULENCE ?= 1
WITH_E4DEBUG ?= 0

# For the ldc2 compiler, we can turn on some extra checking at run-time.
# There will be a computational cost associated with these checks.
WITH_THREAD_SANITIZER ?= 0
WITH_ADDRESS_SANITIZER ?= 0

TECPLOT_BIN_DIR ?= unavailable
TECPLOT_BIN_VERSION_STR = tecplot_unavailable
ifneq ($(TECPLOT_BIN_DIR), unavailable)
    TECPLOT_BIN_VERSION_STR = with_tecplot_binary
endif
TECPLOT_FILES :=
ifneq ($(TECPLOT_BIN_DIR), unavailable)
    TECPLOT_FILES += tecio.d
endif

OPENMPI_DIR := ../extern/OpenMPI
OPENMPI_FILES := $(OPENMPI_DIR)/source/mpi/package.d

MPICH_DIR := ../extern/cray-mpich
MPICH_FILES := $(MPICH_DIR)/mpi.d

PROGRAMS := lmr
SUB_PROGRAMS := lmr-run-steady
SHARE_FILES :=
ETC_FILES := lmr.cfg
PY_PROG_DIR := python-programs
PY_PROGRAMS := lmr-verify

ifeq ($(WITH_MPI),1)
    PARALLEL_FLAVOUR = openmpi
    SUB_PROGRAMS += lmr-mpi-run-steady
    MPI_LIBRARY_DIRS = $(shell mpicc --showme:libdirs)
    MPI_LIB_DIRS_SEARCH = $(foreach d, $(MPI_LIBRARY_DIRS), -L-L$d)
    MPI_DIR = $(OPENMPI_DIR)
    MPI_FILES = $(OPENMPI_FILES)
endif

include eilmer4-files.mk
include lmr-files.mk

EILMER4_FILES = $(LMR4_CORE_FILES) \
	$(LMR4_LUA_FILES) \
	$(LMR4_BC_FILES) \
	$(LMR4_SOLID_FILES) \
	$(LMR4_EFIELD_FILES) \
	$(LMR4_EXTRA_FILES)

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

NML_DIR := ../lib
include $(NML_DIR)/nml_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk
LIBGASF := $(GAS_DIR)/libgasf.a

GRID_DIR := ../grid_utils
include $(GRID_DIR)/grid_utils_files.mk

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

GEOM_DIR := ../geom
include $(GEOM_DIR)/geom_files.mk

GASDYN_DIR := ../gasdyn
include $(GASDYN_DIR)/gasdyn_files.mk

CEQ_DIR := ../extern/ceq/source
LIBCEQ := $(CEQ_DIR)/libceq.a
include $(CEQ_DIR)/ceq_files.mk

NTYPES_DIR := ../ntypes
NTYPES_FILES := $(NTYPES_DIR)/complex.d

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/install/lib

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
BUILD_DIR := ../../build
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
    ifneq ($(TECPLOT_BIN_DIR), unavailable)
        DLINKFLAGS += -L$(TECPLOT_BIN_DIR)/libtecio.so -L$(TECPLOT_BIN_DIR)/sys/libstdc++.so.6
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
    FAST_DFLAGS := -w -g -O2 -release -enable-inlining -boundscheck=off --d-version=flavour_fast -flto=full
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
    ifneq ($(TECPLOT_BIN_DIR), unavailable)
        DLINKFLAGS += -L$(TECPLOT_BIN_DIR)/libtecio.so -L$(TECPLOT_BIN_DIR)/sys/libstdc++.so.6
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
endif

# DIP1008 allows throwing of exceptions in @nogc code. Appeared in 2.079.0, 2018-03-01.
# This rules out the use of gdc for compiling the code.
DFLAGS += -dip1008
DFLAGS += -I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GRID_DIR) -I$(GZIP_DIR)

ifeq ($(DEBUG_CHEM),1)
    DFLAGS += $(DVERSION)debug_chem
endif

ifeq ($(WITH_FPE),1)
    DFLAGS += $(DVERSION)enable_fp_exceptions
endif

ifeq ($(WITH_DVODE),1)
    DFLAGS += $(DVERSION)with_dvode
    DLINKFLAGS += -L-lgfortran
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

default: lmr lmr-run-steady

install: $(PROGRAMS) $(SUB_PROGRAMS) $(PY_PROGRAMS)
	cp $(PROGRAMS) $(INSTALL_DIR)/bin
	cp $(SUB_PROGRAMS) $(INSTALL_DIR)/bin
	cp $(PY_PROGRAMS) $(INSTALL_DIR)/bin
	cp lua-modules/*.lua $(INSTALL_DIR)/lib/
	- mkdir -p $(INSTALL_DIR)/etc
	cp lmr.cfg $(INSTALL_DIR)/etc/

clean:
	- rm $(PROGRAMS)
	- rm $(SUB_PROGRAMS)
	- rm $(PY_PROGRAMS)
	- rm lmr-mpi-run-steady
	- rm *.o
	- rm lmrconfig_with_str_subst.d
	- rm runsteady_shared.d runsteady_openmpi.d
	- cd $(OPENMPI_DIR); make clean
	- cd $(MPICH_DIR); make clean
	- cd $(LUA_DIR); make clean
	- cd $(GEOM_DIR); make clean
	- cd $(GZIP_DIR); make clean
	- cd $(GAS_DIR); make clean; rm -f libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(GRID_DIR); make clean
	- cd $(CEQ_DIR); make clean

$(LIBGASF):
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) libgasf.a

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

$(LIBCEQ):
	cd $(CEQ_DIR); make

lmrconfig_with_str_subst.d: lmrconfig.d
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_FULL_REVISION_STRING_HERE/$(FULL_REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		$< > $@

lmr-run-steady: $(EILMER4_FILES) lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/shared/' $(LMR_CMD)/runsteady.d > runsteady_shared.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ $(DVERSION)runsteady_main \
		$(DVERSION)complex_numbers \
		$(DVERSION)newton_krylov \
		runsteady_shared.d \
		lmrconfig_with_str_subst.d \
		$(EILMER4_FILES) $(LMR_CORE_FILES) $(LMR_LUA_FILES) $(LMR_CMD)/command.d $(LMR_CMD)/cmdhelper.d \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS)

lmr-mpi-run-steady: $(EILMER4_FILES) lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/openmpi/' $(LMR_CMD)/runsteady.d > runsteady_openmpi.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ $(DVERSION)runsteady_main \
		$(DVERSION)mpi_parallel \
		$(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)complex_numbers \
		$(DVERSION)newton_krylov \
		runsteady_openmpi.d \
		lmrconfig_with_str_subst.d \
		$(EILMER4_FILES) $(LMR_CORE_FILES) $(LMR_LUA_FILES) $(LMR_CMD)/command.d $(LMR_CMD)/cmdhelper.d \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS) -L-lmpi


lmr: main.d $(EILMER4_FILES) lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)complex_numbers \
		$(DVERSION)newton_krylov \
		main.d \
		lmrconfig_with_str_subst.d \
		$(EILMER4_FILES) $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS)

#### Python programs
lmr-verify: $(PY_PROG_DIR)/lmr_verify.py
	cp $< $@
	chmod +x $@
