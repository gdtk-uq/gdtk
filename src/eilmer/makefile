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
WITH_CHT ?= 0
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

PROGRAMS := e4shared e4monitor gdtk-module
SHARE_FILES :=

ifeq ($(WITH_COMPLEX_NUMBERS),1)
    PROGRAMS += e4zshared
endif

ifeq ($(WITH_DIAGNOSTICS), 1)
	PROGRAMS += e4shared_diagnostic
endif

ifeq ($(WITH_NK),1)
    PROGRAMS += e4-nk-shared
    PROGRAMS += e4-nk-shared-real
endif

ifeq ($(WITH_CHT),1)
    PROGRAMS += e4-cht-shared
    PROGRAMS += e4-cht-shared-real
endif

ifeq ($(WITH_PIR),1)
    PROGRAMS += e4-pir-shared
endif

ifeq ($(WITH_SSC),1)
    # 2021-05-13 PJ Remove these programs from the list because Kyle is refactoring them.
    # Note that we still need to have the WITH_SSC=1 to build the steady-state solver.
    # PROGRAMS += e4ssc
    # PROGRAMS += e4opt
endif

ifeq ($(WITH_MPI),1)
    PROGRAMS += e4mpi e4loadbalance
    ifeq ($(WITH_COMPLEX_NUMBERS), 1)
        PROGRAMS += e4zmpi
    endif
    ifeq ($(WITH_NK),1)
        PROGRAMS += e4-nk-dist
        PROGRAMS += e4-nk-dist-real
    endif
    ifeq ($(WITH_CHT),1)
        PROGRAMS += e4-cht-dist
        PROGRAMS += e4-cht-dist-real
    endif
    ifeq ($(WITH_PIR),1)
        PROGRAMS += e4-pir-dist
    endif
    MPI_LIBRARY_DIRS = $(shell mpicc --showme:libdirs)
    MPI_LIB_DIRS_SEARCH = $(foreach d, $(MPI_LIBRARY_DIRS), -L-L$d)
    MPI_DIR = $(OPENMPI_DIR)
    MPI_FILES = $(OPENMPI_FILES)
endif

ifeq ($(WITH_MPICH),1)
    PROGRAMS += e4mpi e4loadbalance
    ifeq ($(WITH_COMPLEX_NUMBERS), 1)
        PROGRAMS += e4zmpi
    endif
    ifeq ($(WITH_NK),1)
        PROGRAMS += e4-nk-dist
        PROGRAMS += e4-nk-dist-real
    endif
    ifeq ($(WITH_CHT),1)
        PROGRAMS += e4-cht-dist
        PROGRAMS += e4-cht-dist-real
    endif
    ifeq ($(WITH_PIR),1)
        PROGRAMS += e4-pir-dist
    endif
    MPI_LIBRARY_DIRS = $(shell mpicc -link_info | cut -d' ' -f 3)
    MPI_LIB_DIRS_SEARCH = $(foreach d, $(MPI_LIBRARY_DIRS), -L$d)
    MPI_DIR = $(MPICH_DIR)
    MPI_FILES = $(MPICH_FILES)
endif

ifeq ($(WITH_OPENCL_GPU_CHEM),1)
    PROGRAMS += eilmer-opencl-gpu-chem prep-opencl-gpu-chem-kernel.py
    SHARE_FILES += alpha_qss_kernel_opencl_tmplt.cl
endif

ifeq ($(WITH_CUDA_GPU_CHEM),1)
    PROGRAMS += eilmer-cuda-gpu-chem prep-cuda-gpu-chem-kernel.py
    SHARE_FILES += alpha_qss_kernel_cuda_tmplt.cu libcudakernel.so
endif

ifeq ($(WITH_E4DEBUG),1)
    PROGRAMS += e4shared-debug
endif

EILMER_CORE_FILES := simcore.d simcore_gasdynamic_step.d simcore_solid_step.d \
	simcore_exchange.d simcore_io.d \
	fileutil.d fvcell.d fvvertex.d fvinterface.d \
	fluxcalc.d onedinterp.d flowgradients.d lsqinterp.d \
	conservedquantities.d flowstate.d \
	globalconfig.d globaldata.d \
	block.d fluidblock.d sfluidblock.d ufluidblock.d \
	fluidblockarray.d fluidblockio_old.d fluidblockio_new.d \
	gas_solid_interface.d \
	grid_motion.d grid_motion_udf.d grid_motion_shock_fitting.d \
	history.d loads.d special_block_init.d \
	mass_diffusion.d shockdetectors.d turbulence.d jacobian.d \
	fluidblockio.d celldata.d

EILMER_LUA_FILES := lua_helper.d \
	luaflowsolution.d luaflowstate.d user_defined_source_terms.d

EILMER_EXTRA_FILES := postprocess.d flowsolution.d \
	vtk_writer.d tecplot_writer_classic.d tecplot_writer.d

EILMER_FILES := $(EILMER_CORE_FILES) $(EILMER_LUA_FILES) $(EILMER_EXTRA_FILES)

BC_FILES := bc/package.d \
	bc/boundary_condition.d \
	bc/ghost_cell_effect/package.d \
	bc/ghost_cell_effect/ghost_cell.d \
	bc/ghost_cell_effect/internal_copy_then_reflect.d \
	bc/ghost_cell_effect/flow_state_copy.d \
	bc/ghost_cell_effect/flow_state_copy_from_profile.d \
	bc/ghost_cell_effect/flow_state_copy_from_history.d \
	bc/ghost_cell_effect/synthesise_flow_state.d \
	bc/ghost_cell_effect/extrapolate_copy.d \
	bc/ghost_cell_effect/from_upwind.d \
	bc/ghost_cell_effect/fixed_p.d \
	bc/ghost_cell_effect/fixed_pt.d \
	bc/ghost_cell_effect/from_stagnation.d \
	bc/ghost_cell_effect/full_face_copy.d \
	bc/ghost_cell_effect/mapped_cell_copy.d \
	bc/ghost_cell_effect/gas_solid_full_face_copy.d \
	bc/user_defined_effects.d \
	bc/boundary_flux_effect.d \
	bc/boundary_cell_effect.d \
	bc/boundary_interface_effect.d

SOLID_FILES := solid/solidbc.d \
	solid/solidblock.d \
	solid/solid_ghost_cell.d \
	solid/solid_boundary_flux_effect.d \
	solid/solid_boundary_interface_effect.d \
	solid/solid_gas_full_face_copy.d \
	solid/solid_full_face_copy.d \
	solid/ssolidblock.d \
	solid/solidfvcell.d \
	solid/solidfvinterface.d \
	solid/solidfvvertex.d \
	solid/solidprops.d \
	solid/solidsolution.d \
	solid/solid_udf_source_terms.d \
	solid/luasolidprops.d

LUA_MODULES := configoptions.lua \
	blk_conn.lua \
	bc.lua \
	gridpro.lua \
	grid.lua \
	gridarray.lua \
	flowstate.lua \
	fluidblock.lua \
	fbarray.lua \
	solidblock.lua \
	mpi.lua \
	history.lua \
	zones.lua \
	output.lua \
	sssoptions.lua \
	prep_check.lua

FIELD_FILES := field/field.d \
	field/fieldgmres.d \
	field/fieldconductivity.d \
	field/fieldexchange.d \
	field/fieldderivatives.d \
	field/fieldbc.d

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

OPENCL_DIR := ../extern/opencl.d/source
OPENCL_FILES := $(OPENCL_DIR)/cl.d

CUDA_LIB_DIR ?= /usr/local/cuda/lib64
CUDA_DIR := ../extern/cuda.d/source/cuda_d
CUDA_FILES := $(CUDA_DIR)/cublas_api.d \
	$(CUDA_DIR)/cublasXt.d \
	$(CUDA_DIR)/cudalibxt.d \
	$(CUDA_DIR)/curand.d \
	$(CUDA_DIR)/cublas.d \
	$(CUDA_DIR)/cuComplex.d \
	$(CUDA_DIR)/cuda_profiler_api.d \
	$(CUDA_DIR)/vector_types.d \
	$(CUDA_DIR)/cublas_v2.d \
	$(CUDA_DIR)/cuda.d \
	$(CUDA_DIR)/cuda_runtime_api.d

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

# Some information that we want to embed into the source code.
BUILD_DATE := $(shell date)

REPO_DIR := $(shell cd ../../; pwd)
ifneq ("$(wildcard $(REPO_DIR)/.git)","")
    # Repository exists
    REVISION_STRING := $(shell git rev-parse --short HEAD)
    REVISION_AGE := $(shell git log -1 --format=%cd --date=relative)
    REVISION_DATE := $(shell git log -1 --format=%cd)
else
    REVISION_STRING := unknown
    REVISION_AGE := unknown
    REVISION_DATE := unknown
endif

DFLAGS :=
DLINKFLAGS :=

ifeq ($(DMD), dmd)
    DEBUG_DFLAGS := -w -g -debug -version=flavour_debug
    PROFILE_DFLAGS := -profile -w -g -O -release -boundscheck=off -version=flavour_profile
    FAST_DFLAGS := -w -g -O -release -boundscheck=off -version=flavour_fast
    OF := -of
    DVERSION := -version=
    # DIP1008 allows throwing of exceptions in @nogc code. Appeared in 2.079.0, 2018-03-01.
    DFLAGS += -dip1008 -preview=in
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
    # DIP1008 allows throwing of exceptions in @nogc code. Appeared in 2.079.0, 2018-03-01.
    DFLAGS += -dip1008 -preview=in
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
    DEBUG_DFLAGS := -w -g -link-defaultlib-debug -d-debug -d-version=flavour_debug
    ifeq ($(PLATFORM), macosx)
	DEBUG_DFLAGS += -L-export_dynamic
    else
	DEBUG_DFLAGS += -L-export-dynamic
    endif
    PROFILE_DFLAGS := -fprofile-generate -g -w -O2 -release -enable-inlining -boundscheck=off -d-version=flavour_profile
    FAST_DFLAGS := -w -g -O2 -release -enable-inlining -boundscheck=off -d-version=flavour_fast -flto=full
    # DIP1008 allows throwing of exceptions in @nogc code. Appeared in 2.079.0, 2018-03-01.
    DFLAGS += -dip1008 -preview=in
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
ifeq ($(DMD), gdc-12)
    DEBUG_DFLAGS := -Wall -Og -g -fversion=flavour_debug -fcheck=bounds -fcheck=assert -fdebug
    PROFILE_DFLAGS := -Wall -O2 -frelease -fno-debug -fversion=flavour_profile
    FAST_DFLAGS := -Wall -O2 -frelease -fno-debug -fversion=flavour_fast
    OF := -o
    DVERSION := -fversion=
    # DIP1008 allows throwing of exceptions in @nogc code. Appeared in 2.079.0, 2018-03-01.
    # Needs at least gdc-12.
    DFLAGS += -fpreview=dip1008 -fpreview=in
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

# The compiler needs to be able to find source files in the following directories.
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

default: $(PROGRAMS) ugrid_partition gdtk-module
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Eilmer4 simulation code built."

install: $(PROGRAMS) $(LUA_MODULES) prep.lua prep-grids.lua prep-flow.lua post.lua \
		prep-gas prep-chem chemkin2eilmer prep-kinetics \
		$(SHARE_FILES) ugrid_partition
	- mkdir -p $(INSTALL_DIR)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share
	cp $(PROGRAMS) $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/bin/* $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/lib/lpeg.so $(BUILD_DIR)/lib
	cp prep.lua $(BUILD_DIR)/bin
	cp prep-grids.lua $(BUILD_DIR)/bin
	cp prep-flow.lua $(BUILD_DIR)/bin
	cp post.lua $(BUILD_DIR)/bin
	cp $(LUA_MODULES) $(BUILD_DIR)/lib
	cp -r $(LUNAJSON)/src/* $(BUILD_DIR)/lib
	cp -r ../lib/* $(BUILD_DIR)/lib
	cp $(NML_LUA_MODULES) $(BUILD_DIR)/lib
	cp e4console.tcl $(BUILD_DIR)/bin
	chmod +x $(BUILD_DIR)/bin/e4console.tcl
	cp e4compact.rb $(BUILD_DIR)/bin/e4compact
	chmod +x $(BUILD_DIR)/bin/e4compact
	cp e4-prep-restart.rb $(BUILD_DIR)/bin/e4-prep-restart
	chmod +x $(BUILD_DIR)/bin/e4-prep-restart
	cp e4-prep-parallel.py $(BUILD_DIR)/bin/e4-prep-parallel
	chmod +x $(BUILD_DIR)/bin/e4-prep-parallel
	cp xtdata.rb $(BUILD_DIR)/bin/xtdata.rb
	chmod +x $(BUILD_DIR)/bin/xtdata.rb
	cp e4forces.py $(BUILD_DIR)/bin/e4forces
	chmod +x $(BUILD_DIR)/bin/e4forces
	cp ../grid_utils/ugrid_partition $(BUILD_DIR)/bin
	cp gdtk-module $(BUILD_DIR)/share
ifneq "$(SHARE_FILES)" ""
	cp $(SHARE_FILES) $(BUILD_DIR)/share
endif
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

tarball: install
	- rm dgd-for-mech4480.tgz
	hg archive -p dgd dgd-for-mech4480.tar
	cd ../geom; make INSTALL_DIR=$(INSTALL_DIR) install
	tar -C $(INSTALL_DIR)/.. -r dgdinst -f dgd-for-mech4480.tar
	gzip dgd-for-mech4480.tar
	mv dgd-for-mech4480.tar.gz dgd-for-mech4480.tgz

sync-tarball: tarball
	rsync -av dgd-for-mech4480.tgz cfcfd:/opt/local/share/httpd/htdocs/eilmer/tarball/
	rsync -av dgd-for-mech4480.tgz mango:dropbox/

clean:
	- rm -f *.o
	- rm -f $(DEMO_PROGRAMS) $(PROGRAMS)
	- rm -rf $(BUILD_DIR)/*
	- rm -f test-grid.txt.gz test-flow.txt.gz trace.log trace.def
	- rm -rf plot/ flow/ grid/
	- rm -f main_with_rev_string.d
	- rm -f main_complex_with_rev_string.d
	- rm -f nk_accelerator_with_string_subst.d
	- rm -f mpi_nk_accelerator_with_string_subst.d
	- rm -f steadystate_solver_with_rev_string.d
	- rm -f steadystate_solver_complex_with_rev_string.d
	- rm -f shape_sensitivity_calc_with_rev_string.d
	- rm -f e4opt optdriver_with_rev_string.d
	- rm -f mpi_main_with_rev_string.d
	- rm -f mpi_main_complex_with_rev_string.d
	- rm -f loadbalance_with_rev_string.d
	- rm -f e4mpi e4zmpi e4loadbalance e4ssc e4zshared
	- rm -f e4-nk-shared e4-nk-shared-real e4-nk-dist e4-nk-dist-real
	- rm -f e4-pir-shared e4-pir-shared-real e4-pir-dist e4-pir-dist-real
	- rm -f e4-cht-shared e4-cht-shared-real e4-cht-dist e4-cht-dist-real
	- rm -f e4shared_diagnostic
	- rm -f e4shared-debug
	- rm -f dgd-for-mech4480.tgz
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

prep-gas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) build-prep-gas

ugrid_partition:
	cd $(GRID_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) ugrid_partition

prep-chem:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-prep-chem

chemkin2eilmer:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-chemkin2eilmer

prep-kinetics:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-prep-kinetics

$(LIBGASF):
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) libgasf.a

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local

$(LIBCEQ):
	cd $(CEQ_DIR); make

gdtk-module:
	sed -e 's+PUT_REVISION_STRING_HERE+$(REVISION_STRING)+' \
		-e 's+PUT_COMPILER_NAME_HERE+$(DMD)+' \
		-e 's+PUT_INSTALL_DIR_HERE+$(INSTALL_DIR)+' \
		-e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
		-e 's+PUT_BUILD_DATE_HERE+$(BUILD_DATE)+' \
		-e 's+PUT_REVISION_AGE_HERE+$(REVISION_AGE)+' \
		-e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
		gdtk-module-template > gdtk-module

e4monitor: e4monitor.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)e4monitor e4monitor.d $(DLINKFLAGS)

e4shared: main.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		main.d > main_with_rev_string.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)e4shared \
		$(DVERSION)$(TECPLOT_BIN_VERSION_STR) \
		main_with_rev_string.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

e4shared-debug: main.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		main.d > main_with_rev_string.d
	$(DMD) $(DEBUG_DFLAGS) $(DFLAGS) $(OF)e4shared-debug \
		$(DVERSION)$(TECPLOT_BIN_VERSION_STR) \
		main_with_rev_string.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

e4shared_diagnostic: main.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(MPL_FILES) $(GPERF_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		main.d > main_with_rev_string.d
	python3 $(MPL_DIR)/../python/prebuild.py $(MPL_DIR)/..
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -J$(MPL_DIR)/../views $(OF)e4shared_diagnostic \
		$(DVERSION)diagnostics \
		main_with_rev_string.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(MPL_FILES) $(GPERF_FILES)\
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -I$(MPL_DIR) -I$(GPERF_DIR) -L-lprofiler

e4zshared: main.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		main.d > main_complex_with_rev_string.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)e4zshared \
		$(DVERSION)complex_numbers \
		$(DVERSION)$(TECPLOT_BIN_VERSION_STR) \
		main_complex_with_rev_string.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

e4mpi: main.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		main.d > mpi_main_with_rev_string.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)e4mpi \
		$(DVERSION)mpi_parallel \
		$(MPI_LIB_DIRS_SEARCH) \
		mpi_main_with_rev_string.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -L-lmpi

e4zmpi: main.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		main.d > mpi_main_complex_with_rev_string.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)e4zmpi \
		$(DVERSION)mpi_parallel \
		$(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)complex_numbers \
		mpi_main_complex_with_rev_string.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -L-lmpi

e4-nk-shared: nk_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		nk_accelerator.d > nk_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)complex_numbers \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		nk_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

e4-nk-shared-real: nk_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		nk_accelerator.d > nk_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		nk_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

e4-nk-dist: nk_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		nk_accelerator.d > mpi_nk_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)mpi_parallel $(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)complex_numbers \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		mpi_nk_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -L-lmpi

e4-nk-dist-real: nk_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		nk_accelerator.d > mpi_nk_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)mpi_parallel $(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		mpi_nk_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -L-lmpi

e4-cht-shared: cht_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		cht_accelerator.d > cht_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)complex_numbers \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		cht_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

e4-cht-shared-real: cht_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		cht_accelerator.d > cht_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		cht_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

e4-cht-dist: cht_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		cht_accelerator.d > mpi_cht_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)mpi_parallel $(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)complex_numbers \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		mpi_cht_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -L-lmpi

e4-cht-dist-real: cht_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		cht_accelerator.d > mpi_cht_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)mpi_parallel $(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		mpi_cht_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -L-lmpi

e4-pir-shared: nk_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		nk_accelerator.d > nk_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)pir \
		$(DVERSION)complex_numbers \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		nk_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

e4-pir-shared-real: nk_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		nk_accelerator.d > nk_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)pir \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		nk_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

e4-pir-dist: nk_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		nk_accelerator.d > mpi_nk_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)mpi_parallel $(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)pir \
		$(DVERSION)complex_numbers \
		$(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		mpi_nk_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -L-lmpi

e4-pir-dist-real: nk_accelerator.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	cd $(MPI_DIR); make
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		-e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
		nk_accelerator.d > mpi_nk_accelerator_with_string_subst.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)$@ \
		$(DVERSION)mpi_parallel $(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)pir \
                $(DVERSION)nk_accelerator \
		$(DVERSION)steady_state \
		mpi_nk_accelerator_with_string_subst.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) $(MPI_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -L-lmpi

e4ssc: shape_sensitivity_calc.d shape_sensitivity_core.d steadystate_core.d solid/solid_loose_coupling_update.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		shape_sensitivity_calc.d > shape_sensitivity_calc_with_rev_string.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)e4ssc \
		$(DVERSION)complex_numbers \
		$(DVERSION)shape_sensitivity \
		$(DVERSION)nk_accelerator \
		shape_sensitivity_calc_with_rev_string.d \
		shape_sensitivity_core.d steadystate_core.d solid/solid_loose_coupling_update.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(GAS_LUA_FILES) $(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

e4opt: optdriver.d /
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		optdriver.d > optdriver_with_rev_string.d
	$(DMD)  $(OF)e4opt optdriver_with_rev_string.d

e4loadbalance: loadbalance.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(GAS_LUA_FILES) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GASDYN_FILES) $(GASDYN_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		loadbalance.d > loadbalance_with_rev_string.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)e4loadbalance \
		loadbalance_with_rev_string.d $(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(GAS_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

eilmer-opencl-gpu-chem: main.d opencl_gpu_chem.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) $(OPENCL_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(NM_LUA_FILES) $(GASDYN_FILES) $(GASDYN_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		main.d > main_with_rev_string.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)eilmer-opencl-gpu-chem \
		$(DVERSION)gpu_chem $(DVERSION)opencl_gpu_chem \
		main_with_rev_string.d $(EILMER_FILES) opencl_gpu_chem.d \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) $(OPENCL_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(GAS_LUA_FILES) $(NM_LUA_FILES) $(GASDYN_FILES) $(GASDYN_LUA_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -I$(OPENCL_DIR) -L-L/usr/lib64/nvidia -L-lOpenCL

eilmer-cuda-gpu-chem: main.d cuda_gpu_chem.d alpha_qss_kernel_cuda_tmplt.cu $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) $(CUDA_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(NM_LUA_FILES) $(GASDYN_FILES) $(GASDYN_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		main.d > main_with_rev_string.d
	$(NVCC) --shared -o libcudakernel.so alpha_qss_kernel_cuda_tmplt.cu \
		--compiler-options '-fPIC'
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) $(OF)eilmer-cuda-gpu-chem \
		$(DVERSION)gpu_chem $(DVERSION)cuda_gpu_chem \
		main_with_rev_string.d $(EILMER_FILES) cuda_gpu_chem.d \
		$(BC_FILES) $(SOLID_FILES) $(FIELD_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) $(OPENCL_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(GAS_LUA_FILES) $(NM_LUA_FILES) $(GASDYN_FILES) $(GASDYN_LUA_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS) -I$(CUDA_DIR) $(CUDA_FILES) -L-L$(CUDA_LIB_DIR) \
		-L-lcuda -L-lcudart -L-L. -L-lcudakernel

prep-opencl-gpu-chem-kernel.py: utils/prep-opencl-gpu-chem-kernel.py
	cp utils/prep-opencl-gpu-chem-kernel.py .
	chmod +x prep-opencl-gpu-chem-kernel.py

alpha_qss_kernel_opencl_tmplt.cl: utils/alpha_qss_kernel_opencl_tmplt.cl
	cp utils/alpha_qss_kernel_opencl_tmplt.cl .

prep-cuda-gpu-chem-kernel.py: utils/prep-cuda-gpu-chem-kernel.py
	cp utils/prep-cuda-gpu-chem-kernel.py .
	chmod +x prep-cuda-gpu-chem-kernel.py

alpha_qss_kernel_cuda_tmplt.cu: utils/alpha_qss_kernel_cuda_tmplt.cu
	cp utils/alpha_qss_kernel_cuda_tmplt.cu .

libcudakernel.so: libcudakernel.so
	cp libcudakernel.so ./../lib/
