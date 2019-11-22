# makefile for Eilmer
# Builds main simulation program by default.
# For example:
# make install
#
# As default, we compile with the 'fast' option.
# For the DMD compiler, the 'fast' option requests
# production mode and to include optimisations.
#
# When more detail is required in the error output,
# we recommend compiling with the 'debug' option.
# For the DMD compiler, this invokes all of the assert
# checks. It also does bounds checks on all of the
# array accesses.

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc and ldc2.
DMD ?= dmd
NVCC ?= nvcc
GPP ?= g++

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug

PLATFORM ?= linux
WITH_MPI ?= 0
WITH_MPI_TIMEOUTS ?= 0
WITH_SSS ?= 0
WITH_SSC ?= 0
WITH_OPENCL_GPU_CHEM ?= 0
WITH_CUDA_GPU_CHEM ?= 0
DEBUG_CHEM ?= 0
WITH_COMPLEX_NUMBERS ?= 0
WITH_FPE ?= 0
WITH_DVODE ?= 0
MULTI_SPECIES_GAS ?= 1
MULTI_T_GAS ?= 1
MHD ?= 1
KOMEGA ?= 1

TECPLOT_BIN_DIR ?= unavailable
TECPLOT_BIN_VERSION_STR = tecplot_unavailable
ifneq ($(TECPLOT_BIN_DIR), unavailable)
    TECPLOT_BIN_VERSION_STR = with_tecplot_binary
endif
TECPLOT_FILES :=
ifneq ($(TECPLOT_BIN_DIR), unavailable)
    TECPLOT_FILES += tecio.d
endif



PROGRAMS := e4shared e4monitor
SHARE_FILES :=

ifeq ($(WITH_SSS),1)
    PROGRAMS += e4sss
endif

ifeq ($(WITH_SSC),1)
    PROGRAMS += e4ssc
    PROGRAMS += e4opt
endif

ifeq ($(WITH_MPI),1)
    PROGRAMS += e4mpi e4loadbalance
    MPI_LIBRARY_DIRS = $(shell mpicc --showme:libdirs)
    MPI_LIB_DIRS_SEARCH = $(foreach d, $(MPI_LIBRARY_DIRS), -L-L$d)
endif

ifeq ($(WITH_OPENCL_GPU_CHEM),1)
    PROGRAMS += eilmer-opencl-gpu-chem prep-opencl-gpu-chem-kernel.py
    SHARE_FILES += alpha_qss_kernel_opencl_tmplt.cl
endif

ifeq ($(WITH_CUDA_GPU_CHEM),1)
    PROGRAMS += eilmer-cuda-gpu-chem prep-cuda-gpu-chem-kernel.py
    SHARE_FILES += alpha_qss_kernel_cuda_tmplt.cu libcudakernel.so
endif

ifeq ($(WITH_COMPLEX_NUMBERS),1)
    PROGRAMS += e4zshared
    PROGRAMS += e4zsss
endif

EILMER_CORE_FILES := simcore.d \
	fileutil.d json_helper.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d \
	fluxcalc.d onedinterp.d flowgradients.d lsqinterp.d limiters.d \
	conservedquantities.d flowstate.d \
	globalconfig.d globaldata.d \
	fluidblock.d sfluidblock.d ufluidblock.d \
	gas_solid_interface.d \
	flowsolution.d grid_motion.d \
	shock_fitting_moving_grid.d history.d loads.d special_block_init.d \
	mass_diffusion.d tecplot_writer.d

EILMER_LUA_FILES := luaglobalconfig.d lua_helper.d \
	luaflowsolution.d luaflowstate.d user_defined_source_terms.d

EILMER_EXTRA_FILES := postprocess.d

EILMER_FILES := $(EILMER_CORE_FILES) $(EILMER_LUA_FILES) $(EILMER_EXTRA_FILES)

BC_FILES := bc/package.d \
	bc/boundary_condition.d \
	bc/ghost_cell_effect/package.d \
	bc/ghost_cell_effect/ghost_cell.d \
	bc/ghost_cell_effect/internal_copy_then_reflect.d \
	bc/ghost_cell_effect/flow_state_copy.d \
	bc/ghost_cell_effect/flow_state_copy_from_profile.d \
	bc/ghost_cell_effect/extrapolate_copy.d \
	bc/ghost_cell_effect/fixed_p.d \
	bc/ghost_cell_effect/fixed_pt.d \
	bc/ghost_cell_effect/from_stagnation.d \
	bc/ghost_cell_effect/full_face_copy.d \
	bc/ghost_cell_effect/mapped_cell_copy.d \
	bc/user_defined_effects.d \
	bc/boundary_flux_effect.d \
	bc/boundary_cell_effect.d \
	bc/boundary_interface_effect.d

SOLID_FILES := solid/solidbc.d \
	solid/solidblock.d \
	solid/solid_boundary_flux_effect.d \
	solid/solid_boundary_interface_effect.d \
	solid/solid_loose_coupling_update.d \
	solid/ssolidblock.d \
	solid/solidfvcell.d \
	solid/solidfvinterface.d \
	solid/solidfvvertex.d \
	solid/solidprops.d \
	solid/solidsolution.d \
	solid/solid_udf_source_terms.d \
	solid/luasolidprops.d

LUA_MODULES := blk_conn.lua \
	bc.lua \
	gridpro.lua


UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk
LUA_MODULES += $(UTIL_DIR)/lua_helper.lua \
	$(UTIL_DIR)/simplecsv.lua \
	$(UTIL_DIR)/billig_patch.lua \
	$(UTIL_DIR)/tools.lua

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

NML_DIR := ../lib
include $(NML_DIR)/nml_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk

GRID_DIR := ../grid_utils
include $(GRID_DIR)/grid_utils_files.mk

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

GEOM_DIR := ../geom
include $(GEOM_DIR)/geom_files.mk

GASDYN_DIR := ../gasdyn
include $(GASDYN_DIR)/gasdyn_files.mk

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

OPENMPI_DIR := ../extern/OpenMPI/
OPENMPI_FILES := $(OPENMPI_DIR)/source/mpi/package.d $(OPENMPI_DIR)/source/mpi/util.d

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := $(LUA_DIR)/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib

# The install destination.
INSTALL_DIR ?= $(HOME)/dgdinst

# The build destination sits locally for present
BUILD_DIR := ../../build

REVISION_STRING := $(shell hg identify --id --num --branch --tags)

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
    DLINKFLAGS := -L-L$(LIBLUAPATH) 
    ifeq ($(findstring with_libplot,$(LIBPLOT_VERSION_STR)), with_libplot)
        DFLAGS := $(DFLAGS) $(DVERSION)$(LIBPLOT_VERSION_STR)
        DLINKFLAGS := $(DLINKFLAGS) -L-lplot
    endif
    DLINKFLAGS := $(DLINKFLAGS) -L-llua -L-ldl
    ifneq ($(TECPLOT_BIN_DIR), unavailable)
        DLINKFLAGS += -L$(TECPLOT_BIN_DIR)/libtecio.so -L$(TECPLOT_BIN_DIR)/sys/libstdc++.so.6
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
    OF := -of
    DVERSION := -version=
    DLINKFLAGS := -L-L$(LIBLUAPATH) 
    ifeq ($(findstring with_libplot,$(LIBPLOT_VERSION_STR)), with_libplot)
        DFLAGS := $(DFLAGS) $(DVERSION)$(LIBPLOT_VERSION_STR)
        DLINKFLAGS := $(DLINKFLAGS) -L-lplot
    endif
    DLINKFLAGS := $(DLINKFLAGS) -L-llua -L-ldl
    ifneq ($(TECPLOT_BIN_DIR), unavailable)
        DLINKFLAGS += -L$(TECPLOT_BIN_DIR)/libtecio.so -L$(TECPLOT_BIN_DIR)/sys/libstdc++.so.6
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
        DFLAGS := -w -g -O -release -enable-inlining -boundscheck=off -d-version=flavour_fast -ffast-math -flto=full
    endif
    OF := -of=
    DVERSION := -d-version=
    DLINKFLAGS := -L-L$(LIBLUAPATH)
    #ifeq ($(FLAVOUR), profile)
    #    DLINKFLAGS := $(DLINKFLAGS) -Wl,-fprofile-generate
    #endif
    ifeq ($(findstring with_libplot,$(LIBPLOT_VERSION_STR)), with_libplot)
        DFLAGS := $(DFLAGS) $(DVERSION)$(LIBPLOT_VERSION_STR)
        DLINKFLAGS := $(DLINKFLAGS) -L-lplot
    endif
    DLINKFLAGS := $(DLINKFLAGS) -L-llua -L-ldl
    ifneq ($(TECPLOT_BIN_DIR), unavailable)
        DLINKFLAGS += -L$(TECPLOT_BIN_DIR)/libtecio.so -L$(TECPLOT_BIN_DIR)/sys/libstdc++.so.6
    endif
endif
ifeq ($(DMD), gdc)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -Wall -Og -g -version=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        DFLAGS := -Wall -O2 -frelease -fno-debug -version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -Wall -O2 -frelease -fno-debug -version=flavour_fast
    endif
    OF := -o
    DVERSION := -version=
    DLINKFLAGS := -L$(LIBLUAPATH) 
    ifeq ($(findstring with_libplot,$(LIBPLOT_VERSION_STR)), with_libplot)
        DFLAGS := $(DFLAGS) $(DVERSION)$(LIBPLOT_VERSION_STR)
        DLINKFLAGS := $(DLINKFLAGS) -lplot
    endif
    DLINKFLAGS := $(DLINKFLAGS) $(LIBLUA) -ldl -lcurl
    ifneq ($(TECPLOT_BIN_DIR), unavailable)
        DLINKFLAGS += $(TECPLOT_BIN_DIR)/libtecio.so $(TECPLOT_BIN_DIR)/sys/libstdc++.so.6
    endif
endif
# DIP1008 allows throwing of exceptions in @nogc code.
DFLAGS += -dip1008
DFLAGS += -I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GRID_DIR) -I$(GZIP_DIR)

ifeq ($(WITH_MPI_TIMEOUTS),1)
    DFLAGS += $(DVERSION)mpi_timeouts
endif

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

ifeq ($(KOMEGA),1)
    DFLAGS += $(DVERSION)komega
endif

default: $(PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Eilmer4 simulation code built."

install: $(PROGRAMS) $(LUA_MODULES) prep.lua prep-gas prep-chem chemkin2eilmer \
		$(SHARE_FILES) ugrid_partition
	- mkdir -p $(INSTALL_DIR)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share
	cp $(PROGRAMS) $(BUILD_DIR)/bin
	cp $(LUA_DIR)/bin/* $(BUILD_DIR)/bin
	cp prep.lua $(BUILD_DIR)/bin
	cp post.lua $(BUILD_DIR)/bin
	cp $(LUA_MODULES) $(BUILD_DIR)/lib
	cp $(NML_LUA_MODULES) $(BUILD_DIR)/lib
	cp $(GASDYN_LUA_MODULES) $(BUILD_DIR)/lib
	cp e4console.tcl $(BUILD_DIR)/bin
	chmod +x $(BUILD_DIR)/bin/e4console.tcl
	cp e4compact.rb $(BUILD_DIR)/bin/e4compact
	chmod +x $(BUILD_DIR)/bin/e4compact
	cp xtdata.rb $(BUILD_DIR)/bin/xtdata.rb
	chmod +x $(BUILD_DIR)/bin/xtdata.rb
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
	- rm *.o
	- rm $(DEMO_PROGRAMS) $(PROGRAMS)
	- rm -r $(BUILD_DIR)/*
	- rm test-grid.txt.gz test-flow.txt.gz trace.log trace.def
	- rm -r plot/ flow/ grid/
	- rm main_with_rev_string.d
	- rm main_complex_with_rev_string.d
	- rm steadystate_solver_with_rev_string.d
	- rm steadystate_solver_complex_with_rev_string.d
	- rm shape_sensitivity_calc_with_rev_string.d
	- rm e4opt optdriver_with_rev_string.d
	- rm mpi_main_with_rev_string.d
	- rm loadbalance_with_rev_string.d
	- rm e4mpi e4loadbalance e4ssc e4sss e4zshared e4zsss
	- rm dgd-for-mech4480.tgz
	- cd $(OPENMPI_DIR); make clean
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(GRID_DIR); make clean

prep-gas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) build-prep-gas

ugrid_partition:
	cd $(GRID_DIR); make BUILD_DIR=$(BUILD_DIR) ugrid_partition


prep-chem:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-prep-chem

chemkin2eilmer:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-chemkin2eilmer


$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local

e4monitor: e4monitor.d
	$(DMD) $(DFLAGS) $(OF)e4monitor e4monitor.d $(DLINKFLAGS)

e4shared: main.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		main.d > main_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4shared \
		$(DVERSION)$(TECPLOT_BIN_VERSION_STR) \
		main_with_rev_string.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(DLINKFLAGS)

e4zshared: main.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		main.d > main_complex_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4zshared \
		$(DVERSION)complex_numbers \
		$(DVERSION)$(TECPLOT_BIN_VERSION_STR) \
		main_complex_with_rev_string.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(DLINKFLAGS)

e4mpi: main.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(NM_LUA_FILES)
	cd $(OPENMPI_DIR); gen/setup.sh
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		main.d > mpi_main_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4mpi \
		$(DVERSION)mpi_parallel \
		-I$(OPENMPI_DIR) $(MPI_LIB_DIRS_SEARCH) \
		mpi_main_with_rev_string.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(NM_LUA_FILES) $(OPENMPI_FILES) \
		$(DLINKFLAGS) -L-lmpi

e4sss: steadystate_solver.d steadystate_core.d shape_sensitivity_core.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(GASDYN_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		steadystate_solver.d > steadystate_solver_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4sss \
		$(DVERSION)steady_state \
		$(DVERSION)shape_sensitivity \
		steadystate_solver_with_rev_string.d steadystate_core.d shape_sensitivity_core.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(GASDYN_FILES) $(NM_LUA_FILES) \
		$(DLINKFLAGS)

e4zsss: steadystate_solver.d steadystate_core.d shape_sensitivity_core.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		steadystate_solver.d > steadystate_solver_complex_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4zsss \
		$(DVERSION)complex_numbers \
		$(DVERSION)steady_state \
		$(DVERSION)shape_sensitivity \
		steadystate_solver_complex_with_rev_string.d steadystate_core.d shape_sensitivity_core.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(NM_LUA_FILES) $(TECPLOT_FILES) \
		$(DLINKFLAGS)

e4ssc: shape_sensitivity_calc.d shape_sensitivity_core.d steadystate_core.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(GASDYN_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		shape_sensitivity_calc.d > shape_sensitivity_calc_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4ssc \
		$(DVERSION)complex_numbers \
		$(DVERSION)shape_sensitivity \
		$(DVERSION)steady_state \
		shape_sensitivity_calc_with_rev_string.d \
		shape_sensitivity_core.d steadystate_core.d \
		$(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(GASDYN_FILES) $(NM_LUA_FILES) \
		$(DLINKFLAGS)

e4opt: optdriver.d /
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		optdriver.d > optdriver_with_rev_string.d
	$(DMD)  $(OF)e4opt optdriver_with_rev_string.d

e4loadbalance: loadbalance.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GASDYN_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		loadbalance.d > loadbalance_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4loadbalance \
		loadbalance_with_rev_string.d $(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GASDYN_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(NM_LUA_FILES) \
		$(DLINKFLAGS)

eilmer-opencl-gpu-chem: main.d opencl_gpu_chem.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(LIBLUA) $(GZIP_FILES) $(OPENCL_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(NM_LUA_FILES) $(GASDYN_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		main.d > main_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)eilmer-opencl-gpu-chem \
		$(DVERSION)gpu_chem $(DVERSION)opencl_gpu_chem \
		main_with_rev_string.d $(EILMER_FILES) opencl_gpu_chem.d \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(GZIP_FILES) $(OPENCL_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(NM_LUA_FILES) $(GASDYN_FILES) \
		$(DLINKFLAGS) -I$(OPENCL_DIR) -L-L/usr/lib64/nvidia -L-lOpenCL

eilmer-cuda-gpu-chem: main.d cuda_gpu_chem.d alpha_qss_kernel_cuda_tmplt.cu $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(LIBLUA) $(GZIP_FILES) $(CUDA_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES) $(NM_LUA_FILES) $(GASDYN_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		main.d > main_with_rev_string.d
	$(NVCC) --shared -o libcudakernel.so alpha_qss_kernel_cuda_tmplt.cu \
		--compiler-options '-fPIC'
	$(DMD) $(DFLAGS) $(OF)eilmer-cuda-gpu-chem \
		$(DVERSION)gpu_chem $(DVERSION)cuda_gpu_chem \
		main_with_rev_string.d $(EILMER_FILES) cuda_gpu_chem.d \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(GZIP_FILES) $(OPENCL_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(NM_LUA_FILES) $(GASDYN_FILES) \
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

#--------------------------------------------------------------------------
# Fortran bits -- we need a bit of a hack

$(GAS_DIR)/ideal_air_fortran.o: $(GAS_DIR)/ideal_air_fortran.f
	cd $(GAS_DIR); gfortran -c -ffree-form -fPIC ideal_air_fortran.f

$(KINETICS_DIR)/pseudosp_solve_ode.o: $(KINETICS_DIR)/pseudosp_solve_ode.f90
	cd $(KINETICS_DIR); gfortran -c pseudosp_solve_ode.f90

$(KINETICS_DIR)/pseudosp_rhs_ode.o: $(KINETICS_DIR)/pseudosp_rhs_ode.f90
	cd $(KINETICS_DIR); gfortran -c pseudosp_rhs_ode.f90

$(KINETICS_DIR)/dvode.o: $(KINETICS_DIR)/dvode.f
	cd $(KINETICS_DIR); gfortran -c dvode.f
