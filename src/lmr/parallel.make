DMD ?= ldc2
NVCC ?= nvcc
GPP ?= g++
OPT_OUT_MPI ?= 0
MPI_IMPLEMENTATION := OpenMPI
WITH_MPICH ?= 0

FLAVOUR ?= debug
WITH_FPE ?= 1

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
WITH_CHECK_JAC ?= 0

WITH_THREAD_SANITIZER ?= 0
WITH_ADDRESS_SANITIZER ?= 0

OPENMPI_DIR := ../extern/OpenMPI
OPENMPI_FILES := $(OPENMPI_DIR)/source/mpi/package.d
MPI_DIR = $(OPENMPI_DIR)
MPI_FILES = $(OPENMPI_FILES)
MPICH_DIR := ../extern/cray-mpich
MPICH_FILES := $(MPICH_DIR)/mpi.d

PROGRAMS := lmr lmr-debug gdtk-module
SUB_PROGRAMS := lmr-run lmrZ-run
ifeq ($(OPT_OUT_MPI), 0)
	SUB_PROGRAMS += lmr-mpi-run lmrZ-mpi-run
endif
ifeq ($(WITH_MPICH),1)
    MPI_LIBRARY_DIRS = $(shell mpicc -link_info | cut -d' ' -f 3)
    MPI_LIB_DIRS_SEARCH = $(foreach d, $(MPI_LIBRARY_DIRS), -L$$d)
    MPI_DIR = $(MPICH_DIR)
    MPI_FILES = $(MPICH_FILES)
	MPI_IMPLEMENTATION := MPICH
endif
ifeq ($(WITH_CHECK_JAC), 1)
	SUB_PROGRAMS += lmr-check-jacobian lmrZ-check-jacobian
endif

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

DFLAGS += -I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GRID_DIR) -I$(GZIP_DIR) -I$(CEQ_DIR) -I$(DYAML_DIR)/.. -I$(TINYENDIAN_DIR)

INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR := $(abspath ../../build)
OBJ_DIR := $(abspath $(BUILD_DIR)/obj)
BUILD_DATE := $(shell date)
REVISION_STRING := $(shell git rev-parse --short HEAD)
FULL_REVISION_STRING := $(shell git rev-parse HEAD)
REVISION_AGE := $(shell git log -1 --format=%cd --date=relative)
REVISION_DATE := $(shell git log -1 --format=%cd)
REPO_DIR := $(shell cd ../../; pwd)

CEQ_DIR := ../extern/ceq/source
CEQ_BUILD_DIR := $(abspath $(OBJ_DIR)/extern/ceq/source)
LIBCEQ := $(CEQ_BUILD_DIR)/libceq.a
include $(CEQ_DIR)/ceq_files_parallel.mk

LUA_DIR := $(abspath ../../extern/lua-5.4.3/src)
LUA_BUILD_DIR := $(abspath $(BUILD_DIR)/extern/lua)
LIBLUA := $(LUA_BUILD_DIR)/liblua.a
LIBLUAPATH := $(LUA_BUILD_DIR)

ifeq ($(DMD), ldc2)
    DEBUG_DFLAGS := -w -g --d-debug --d-version=flavour_debug
    PROFILE_DFLAGS := -fprofile-generate -g -w -O2 -release -enable-inlining -boundscheck=off --d-version=flavour_profile
    ifeq ($(PLATFORM), macosx)
        FAST_DFLAGS := -w -g -O1 -release -enable-inlining -boundscheck=off --d-version=flavour_fast
    else
        FAST_DFLAGS := -w -g -O2 -release -enable-inlining -boundscheck=off --d-version=flavour_fast
    endif
    ifeq ($(PLATFORM), linux)
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

    # pass LINKER_FLAGS="--linker=mold" for fast linking
    DLINKFLAGS := $(LINKER_FLAGS)
    ifeq ($(WITH_THREAD_SANITIZER), 1)
        DLINKFLAGS := $(DLINKFLAGS) -fsanitize=thread
    endif
    ifeq ($(WITH_ADDRESS_SANITIZER), 1)
        DLINKFLAGS := $(DLINKFLAGS) -fsanitize=address
    endif
    DLINKFLAGS := $(DLINKFLAGS) -L-ldl
    ifeq ($(PLATFORM), linux)
		DLINKFLAGS += -L-Wl,-E
	endif
    ifeq ($(PLATFORM), macosx)
		DLINKFLAGS += -L-ld_classic
    endif
endif

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

ifeq ($(PLATFORM), macosx)
	DFLAGS += $(DVERSION)macosx
endif

DFLAGS += -dip1008 -I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GRID_DIR) -I$(GZIP_DIR)
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

DFLAGS += $(DVERSION)newton_krylov

SRC_DIR = $(abspath ../)
CONFIG_BUILD_DIR = $(BUILD_DIR)/config

LMR_SRCS = $(abspath main.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
           $(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
           $(GEOM_FILES) $(GRID_FILES) \
           $(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
           $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
           $(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
           $(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
           $(DYAML_FILES) $(TINYENDIAN_FILES))

LMR_OBJS = $(patsubst $(SRC_DIR)/%.d,$(OBJ_DIR)/%.o,$(LMR_SRCS))

$(LMR_OBJS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.d | $(LIBLUA) $(LIBCEQ)
	@mkdir -p $(dir $@)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -c -of=$@ $<

$(CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d: lmrconfig.d
	@mkdir -p $(dir $@)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
	    -e 's/PUT_FULL_REVISION_STRING_HERE/$(FULL_REVISION_STRING)/' \
	    -e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
	    -e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
	    -e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
	    $< > $@

$(CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.o: $(CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d
	@mkdir -p $(dir $@)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -c -of=$@ $<

# $(CONFIG_BUILD_DIR)/runsim_shared_run.d: $(LMR_CMD)/runsim.d $(CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d
# 	@mkdir -p $(dir $@)
# 	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/shared/' \
# 	    -e 's/PUT_NUMBER_TYPE_HERE/real/' \
# 	    $< > $@
# 
# $(CONFIG_BUILD_DIR)/runsim_shared_z.d: $(LMR_CMD)/runsim.d $(CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d
# 	@mkdir -p $(dir $@)
# 	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/shared/' \
# 	    -e 's/PUT_NUMBER_TYPE_HERE/complex/' \
# 	    $< > $@
# 
# $(CONFIG_BUILD_DIR)/runsim_mpi_run.d: $(LMR_CMD)/runsim.d $(CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d
# 	@mkdir -p $(dir $@)
# 	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/$(MPI_IMPLEMENTATION)/' \
# 	    -e 's/PUT_NUMBER_TYPE_HERE/real/' \
# 	    $< > $@
# 
# $(CONFIG_BUILD_DIR)/runsim_mpi_z.d: $(LMR_CMD)/runsim.d $(CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d
# 	@mkdir -p $(dir $@)
# 	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/$(MPI_IMPLEMENTATION)/' \
# 	    -e 's/PUT_NUMBER_TYPE_HERE/complex/' \
# 	    $< > $@


.PHONY: lmr
lmr: $(BUILD_DIR)/lmr

$(BUILD_DIR)/lmr: $(LMR_OBJS) $(LIBCEQ) $(LIBLUA) $(CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.o
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -of=$@ $(LMR_OBJS) $(CONFIG_BUILD_DIR)/lmrconfig_with_str_subst $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

lmr-debug: $(LMR_OBJS) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(DEBUG_DFLAGS) $(DFLAGS) -of=$@ $(LMR_OBJS) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

LMR_RUN_SRCS = runsim_shared_run.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
               $(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
               $(GEOM_FILES) $(GRID_FILES) \
               $(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
               $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
               $(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
               $(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
               $(DYAML_FILES) $(TINYENDIAN_FILES)
LMR_RUN_OBJS = $(patsubst %.d,$(OBJ_DIR)/%.o,$(LMR_RUN_SRCS))

lmr-run: runsim_shared_run.d $(LMR_RUN_OBJS) $(LIBCEQ) $(LIBLUA)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -of=$@ -Drun_main $(LMR_RUN_OBJS) $(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)


$(LIBLUA):
	$(MAKE) LUA_BUILD_DIR=$(LUA_BUILD_DIR) -C $(LUA_DIR) -f parallel.make PLATFORM=$(PLATFORM)
	$(MAKE) LP_BUILD_DIR=$(LUA_BUILD_DIR) -C $(LUA_DIR)/../lpeg-1.0.2 -f parallel.make $(PLATFORM)

$(LIBCEQ):
	$(MAKE) CEA_BUILD_DIR=$(CEQ_BUILD_DIR) -C $(CEQ_DIR) -f parallel.make

$(LIBGASF):
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) libgasf.a

lmr-verify: $(PY_PROG_DIR)/lmr_verify.py
	cp $< $@
	chmod +x $@

prep-gas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) PLATFORM=$(PLATFORM) build-prep-gas

ugrid_partition:
	- mkdir -p $(BUILD_DIR)/bin
	cd $(GRID_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) PLATFORM=$(PLATFORM) ugrid_partition
	cd $(GRID_DIR); cp ugrid_partition $(BUILD_DIR)/bin

prep-chem:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) PLATFORM=$(PLATFORM) build-prep-chem

chemkin2eilmer:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) PLATFORM=$(PLATFORM) build-chemkin2eilmer

prep-kinetics:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) PLATFORM=$(PLATFORM) build-prep-kinetics

gdtk-module:
	sed -e 's+PUT_REVISION_STRING_HERE+$(REVISION_STRING)+' \
	    -e 's+PUT_COMPILER_NAME_HERE+$(DMD)+' \
	    -e 's+PUT_INSTALL_DIR_HERE+$(INSTALL_DIR)+' \
	    -e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
	    -e 's+PUT_BUILD_DATE_HERE+$(BUILD_DATE)+' \
	    -e 's+PUT_REVISION_AGE_HERE+$(REVISION_AGE)+' \
	    -e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
	    ../eilmer/gdtk-module-template > gdtk-module

.PHONY: default install clean prep-gas ugrid_partition prep-chem chemkin2eilmer prep-kinetics gdtk-module

default: $(PROGRAMS) $(SUB_PROGRAMS)

# TODO: FIXME
install: $(PROGRAMS) $(SUB_PROGRAMS) $(PY_PROGRAMS) $(AUX_PROGRAMS)
	- mkdir -p $(INSTALL_DIR)/bin $(INSTALL_DIR)/lib $(INSTALL_DIR)/etc $(INSTALL_DIR)/share
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
	- rm -f $(OBJ_DIR)/*.o
	- rm -rf $(BUILD_DIR)/*
	- rm -f $(PROGRAMS) $(SUB_PROGRAMS) $(PY_PROGRAMS)
	- rm -f lmrconfig_with_str_subst.d runsim_shared_run.d runsim_shared_z.d runsim_mpi_run.d runsim_mpi_z.d
	- cd $(OPENMPI_DIR); make clean
	- cd $(MPICH_DIR); make clean
	- cd $(LUA_DIR); make clean
	- cd $(GEOM_DIR); make clean
	- cd $(GZIP_DIR); make clean
	- cd $(GAS_DIR); make clean; rm -f libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(GRID_DIR); make clean
	- cd $(CEQ_DIR); make clean
