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

# PLATFORM options are linux, macosx
PLATFORM ?= linux

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

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/install/lib

LUNAJSON := ../../extern/lunajson

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
endif
ifeq ($(DMD), ldmd2)
    DEBUG_DFLAGS := -w -g -debug -version=flavour_debug
    PROFILE_DFLAGS := -profile -w -O2 -release -inline -boundscheck=off -version=flavour_profile
    FAST_DFLAGS := -w -g -O2 -release -inline -boundscheck=off -version=flavour_fast
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
endif
ifeq ($(DMD), ldc2)
    # -fprofile-generate will result in profraw files being written
    # that may be viewed, showing the top 10 functions with internal block counts
    # llvm-profdata show -text -topn=10 <profraw-file>
    DEBUG_DFLAGS := -w -g -d-debug -d-version=flavour_debug
    PROFILE_DFLAGS := -fprofile-generate -g -w -O2 -release -enable-inlining -boundscheck=off -d-version=flavour_profile
    FAST_DFLAGS := -w -g -O2 -release -enable-inlining -boundscheck=off -d-version=flavour_fast -flto=full
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
endif
ifeq ($(DMD), gdc-8)
    DEBUG_DFLAGS := -Wall -Og -g -fversion=flavour_debug
    PROFILE_DFLAGS := -Wall -O2 -frelease -fno-debug -fversion=flavour_profile
    FAST_DFLAGS := -Wall -O2 -frelease -fno-debug -fversion=flavour_fast
    OF := -o
    DVERSION := -fversion=
    ifeq ($(findstring with_libplot,$(LIBPLOT_VERSION_STR)), with_libplot)
        DFLAGS := $(DFLAGS) $(DVERSION)$(LIBPLOT_VERSION_STR)
        DLINKFLAGS := $(DLINKFLAGS) -lplot
    endif
    DLINKFLAGS := $(DLINKFLAGS) -ldl -lcurl
    ifneq ($(TECPLOT_BIN_DIR), unavailable)
        DLINKFLAGS += $(TECPLOT_BIN_DIR)/libtecio.so $(TECPLOT_BIN_DIR)/sys/libstdc++.so.6
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
# gdc-8 built on 2.069.2 with updates to 2016-01-15
# gdc-9 built on 2.076.0 which appeared 2017-09-01
# gdc-10 will back-port static foreach, which we also use.
# Expect gdc-10 to be out in May 2020.
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

install:
	cp lmr /work/gdtkinst/bin
	cp lmr-run-steady /work/gdtkinst/bin
	cp lua-modules/*.lua /work/gdtkinst/lib/
	- mkdir -p /work/gdtkinst/etc
	cp lmr.cfg /work/gdtkinst/etc/

clean:
	- rm $(PROGRAMS)
	- rm $(SUB_PROGRAMS)
	- rm lmr.o
	- rm lmr_with_rev_string.d
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
	cd $(LUA_DIR); make $(PLATFORM) local

$(LIBCEQ):
	cd $(CEQ_DIR); make

lmr-run-steady: $(EILMER4_FILES) $(LMR_CORE_FILES) $(LMR_CMD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		$(LMR_CMD)/runsteady.d > main_with_rev_string.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ $(DVERSION)runsteady_main \
		$(DVERSION)complex_numbers \
		$(DVERSION)newton_krylov \
		main_with_rev_string.d \
		$(EILMER4_FILES) $(LMR_CORE_FILES) $(LMR_CMD)/command.d \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

lmr: main.d $(EILMER4_FILES) $(LMR_CORE_FILES) $(LMR_CMD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		main.d > main_with_rev_string.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)complex_numbers \
		$(DVERSION)newton_krylov \
		main_with_rev_string.d \
		$(EILMER4_FILES) $(LMR_CORE_FILES) $(LMR_CMD_FILES) \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)


