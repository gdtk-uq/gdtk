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
FLAVOUR ?= debug
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

PROGRAMS := lmr gdtk-module
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
	# For note on "-Wl,-E", see: https://web.tecgraf.puc-rio.br/~lhf/ftp/lua/install.html#adding
	DLINKFLAGS := $(DLINKFLAGS) -L-Wl,-E -L-ldl
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
    FAST_DFLAGS := -w -g -O2 -release -enable-inlining -boundscheck=off --d-version=flavour_fast
    ifeq ($(PLATFORM), linux)
	# Full link-time optimization does not play nicely on macOS
	FAST_DFLAGS += -flto=full
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
DFLAGS += -I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GRID_DIR) -I$(GZIP_DIR)

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

default: $(PROGRAMS) $(SUB_PROGRAMS)

install: $(PROGRAMS) $(SUB_PROGRAMS) $(PY_PROGRAMS) $(AUX_PROGRAMS)
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
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

clean:
	- rm *.o
	- rm -rf $(BUILD_DIR)/*
	- rm $(PROGRAMS)
	- rm $(SUB_PROGRAMS)
	- rm $(PY_PROGRAMS)
	- rm lmrconfig_with_str_subst.d
	- rm runsim_shared.d runsim_mpi.d
	- rm -f gdtk-module
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

lmr: main.d lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)newton_krylov \
		main.d \
		lmrconfig_with_str_subst.d \
		$(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
		$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS)

lmr-run: lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/shared/' \
	    -e 's/PUT_NUMBER_TYPE_HERE/real/' \
	$(LMR_CMD)/runsim.d > runsim_shared.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ $(DVERSION)run_main \
		$(DVERSION)newton_krylov \
		runsim_shared.d \
		lmrconfig_with_str_subst.d \
		$(LMR_CORE_FILES) $(LMR_LUA_FILES) $(LMR_CMD)/command.d $(LMR_CMD)/cmdhelper.d \
		$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
		$(GEOM_FILES) $(GRID_FILES) \
		-i=$(GAS_PACKAGE) -I=$(GAS_DIR) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS)

lmrZ-run: lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/shared/' \
	    -e 's/PUT_NUMBER_TYPE_HERE/complex/' \
	$(LMR_CMD)/runsim.d > runsim_shared.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ $(DVERSION)run_main \
		$(DVERSION)newton_krylov \
		$(DVERSION)complex_numbers \
		runsim_shared.d \
		lmrconfig_with_str_subst.d \
		$(LMR_CORE_FILES) $(LMR_LUA_FILES) $(LMR_CMD)/command.d $(LMR_CMD)/cmdhelper.d \
		$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
		$(GEOM_FILES) $(GRID_FILES) \
		-i=$(GAS_PACKAGE) -I=$(GAS_DIR) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS)

lmr-mpi-run: lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/$(MPI_IMPLEMENTATION)/' \
	    -e 's/PUT_NUMBER_TYPE_HERE/real/' \
	$(LMR_CMD)/runsim.d > runsim_mpi.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ $(DVERSION)run_main \
		$(DVERSION)mpi_parallel \
		$(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)newton_krylov \
		runsim_mpi.d \
		lmrconfig_with_str_subst.d \
		$(LMR_CORE_FILES) $(LMR_LUA_FILES) $(LMR_CMD)/command.d $(LMR_CMD)/cmdhelper.d \
		$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS) -L-lmpi

lmrZ-mpi-run: lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/$(MPI_IMPLEMENTATION)/' \
	    -e 's/PUT_NUMBER_TYPE_HERE/complex/' \
	$(LMR_CMD)/runsim.d > runsim_mpi.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ $(DVERSION)run_main \
		$(DVERSION)mpi_parallel \
		$(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)newton_krylov \
		$(DVERSION)complex_numbers \
		runsim_mpi.d \
		lmrconfig_with_str_subst.d \
		$(LMR_CORE_FILES) $(LMR_LUA_FILES) $(LMR_CMD)/command.d $(LMR_CMD)/cmdhelper.d \
		$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS) -L-lmpi

lmr-check-jacobian: lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)newton_krylov \
		$(LMR_CMD)/checkjacobian.d \
		lmrconfig_with_str_subst.d \
		$(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
		$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS)

lmrZ-check-jacobian: lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)newton_krylov \
		$(DVERSION)complex_numbers \
		$(LMR_CMD)/checkjacobian.d \
		lmrconfig_with_str_subst.d \
		$(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
		$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS)

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

prep-chem:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) PLATFORM=$(PLATFORM) build-prep-chem

chemkin2eilmer:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) PLATFORM=$(PLATFORM) build-chemkin2eilmer

prep-kinetics:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) PLATFORM=$(PLATFORM) build-prep-kinetics

gdtk-module:
	sed -e 's+PUT_REVISION_STRING_HERE+$(REVISION_STRING)+' \
		-e 's+PUT_COMPILER_NAME_HERE+$(DMD)+' \
		-e 's+PUT_INSTALL_DIR_HERE+$(INSTALL_DIR)+' \
		-e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
		-e 's+PUT_BUILD_DATE_HERE+$(BUILD_DATE)+' \
		-e 's+PUT_REVISION_AGE_HERE+$(REVISION_AGE)+' \
		-e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
		../eilmer/gdtk-module-template > gdtk-module
