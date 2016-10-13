# makefile for Eilmer
# Builds main simulation program by default.
# For example:
# make TARGET=fast install

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd

# TARGET options are debug, fast, profile
# Flags for each compiler will be determined on this option.
TARGET ?= debug

PLATFORM ?= linux
WITH_GPU_CHEM ?= 0

PROGRAMS := e4shared e4sss
SHARE_FILES := 

ifeq ($(WITH_GPU_CHEM),1)
    PROGRAMS += eilmer-gpu-chem prep-gpu-chem-kernel.py
    SHARE_FILES += alpha_qss_kernel.tmplt
endif

DEMO_PROGRAMS := globalconfig_demo flowstate_demo fv_demo block_demo \
	luaflowstate_demo luaglobalconfig_demo

EILMER_FILES := simcore.d fileutil.d json_helper.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		fluxcalc.d onedinterp.d flowgradients.d lsqinterp.d \
		conservedquantities.d flowstate.d luaflowstate.d \
		globalconfig.d luaglobalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		lua_helper.d user_defined_source_terms.d gas_solid_interface.d \
		postprocess.d flowsolution.d luaflowsolution.d grid_motion.d \
		block_moving_grid.d history.d

BC_FILES := bc/bc.d \
	bc/ghost_cell_effect.d \
	bc/user_defined_effects.d \
	bc/boundary_flux_effect.d \
	bc/boundary_interface_effect.d

SOLID_FILES := solid/solidbc.d \
	solid/solidblock.d \
	solid/solid_boundary_flux_effect.d \
	solid/solid_boundary_interface_effect.d \
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

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

GEOM_DIR := ../geom
include $(GEOM_DIR)/geom_files.mk

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

OPENCL_DIR := ../extern/opencl.d/source
OPENCL_FILES := $(OPENCL_DIR)/cl.d

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := $(LUA_DIR)/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib

# The install destination.
INSTALL_DIR ?= $(HOME)/dgdinst

# The build destination sits locally for present
BUILD_DIR := ../../build

REVISION_STRING := $(shell hg identify --id --num --branch --tags)

ifeq ($(DMD), dmd)
    ifeq ($(TARGET), debug)
        DFLAGS := -w -gc -debug
    endif
    ifeq ($(TARGET), profile)
        DFLAGS := -profile -vgc -w -O -release -inline -boundscheck=off
    endif
    ifeq ($(TARGET), fast)
        DFLAGS := -w -O -release -inline -boundscheck=off
    endif
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), ldmd2)
    ifeq ($(TARGET), debug)
        DFLAGS := -w -gc -debug
    endif
    ifeq ($(TARGET), profile)
        DFLAGS := -profile -vgc -w -O -release -inline -boundscheck=off
    endif
    ifeq ($(TARGET), fast)
        DFLAGS := -w -O -release -inline -boundscheck=off
    endif
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), gdc)
    ifeq ($(TARGET), debug)
        DFLAGS := -Wall -Og -g
    endif
    ifeq ($(TARGET), profile)
        DFLAGS := -Wall -O2 -frelease -fno-debug
    endif
    ifeq ($(TARGET), fast)
        DFLAGS := -Wall -O2 -frelease -fno-debug
    endif
    OF := -o
    DLINKFLAGS := -L$(LIBLUAPATH) $(LIBLUA) -ldl -lcurl
endif
DFLAGS += -I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GZIP_DIR)

default: $(PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Eilmer4 simulation code built."

install: $(PROGRAMS) $(LUA_MODULES) prep.lua prep-gas prep-chem $(SHARE_FILES)
	- mkdir -p $(INSTALL_DIR)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share
	cp $(PROGRAMS) $(BUILD_DIR)/bin
	cp $(LUA_DIR)/bin/* $(BUILD_DIR)/bin
	cp prep.lua $(BUILD_DIR)/bin
	cp $(LUA_MODULES) $(BUILD_DIR)/lib
ifneq "$(SHARE_FILES)" ""
	cp $(SHARE_FILES) $(BUILD_DIR)/share
endif
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)


demo: $(DEMO_PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Demo codes built."

clean:
	- rm *.o
	- rm $(DEMO_PROGRAMS) $(PROGRAMS)
	- rm -r $(BUILD_DIR)/*
	- rm test-grid.txt.gz test-flow.txt.gz trace.log trace.def
	- rm -r plot/ flow/ grid/
	- rm main_with_rev_string.d steadystate_solver_with_rev_string.d
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a

prep-gas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) build-prep-gas 

prep-chem:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-prep-chem

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local

globalconfig_demo: globalconfig_demo.d globalconfig.d fvcore.d \
	$(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(GZIP_FILES)
	$(DMD) $(DFLAGS) $(OF)globalconfig_demo \
		globalconfig_demo.d globalconfig.d fvcore.d \
		$(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)  

luaglobalconfig_demo: luaglobalconfig_demo.d luaglobalconfig.d globalconfig.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d conservedquantities.d \
	flowstate.d json_helper.d \
	$(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(GZIP_FILES)
	$(DMD) $(DFLAGS) $(OF)luaglobalconfig_demo \
		luaglobalconfig_demo.d luaglobalconfig.d globalconfig.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d conservedquantities.d \
		flowstate.d json_helper.d \
		$(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(DLINKFLAGS)  

flowstate_demo: flowstate_demo.d flowstate.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d conservedquantities.d \
	globalconfig.d json_helper.d \
	$(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(GZIP_FILES)
	$(DMD) $(DFLAGS) $(OF)flowstate_demo flowstate_demo.d \
		fvcore.d flowstate.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d globalconfig.d json_helper.d \
		$(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(DLINKFLAGS)

luaflowstate_demo: luaflowstate_demo.d $(EILMER_FILES) \
	$(GEOM_FILES) $(BC_FILES) $(SOLID_FILES) $(GAS_FILES) \
	$(LIBLUA) $(UTIL_FILES) $(NM_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(GZIP_FILES)
	$(DMD) $(DFLAGS) $(OF)luaflowstate_demo luaflowstate_demo.d \
		$(EILMER_FILES) \
		$(GEOM_FILES) $(BC_FILES) $(SOLID_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(DLINKFLAGS)

fv_demo: fv_demo.d $(EILMER_FILES) $(SOLID_FILES) \
	$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) 
	$(DMD) $(DFLAGS) $(OF)fv_demo fv_demo.d \
		$(EILMER_FILES) $(SOLID_FILES) \
		$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(DLINKFLAGS)

block_demo: block_demo.d $(EILMER_FILES) $(SOLID_FILES) \
	$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) 
	$(DMD) $(DFLAGS) $(OF)block_demo block_demo.d \
		$(EILMER_FILES) $(SOLID_FILES) \
		$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) \
		$(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(DLINKFLAGS)

e4shared: main.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(NM_FILES) $(UTIL_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' main.d > main_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4shared \
		main_with_rev_string.d $(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(DLINKFLAGS)

e4sss: steadystate_solver.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(NM_FILES) $(UTIL_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' steadystate_solver.d > steadystate_solver_with_rev_string.d
	$(DMD) -version=steady_state $(DFLAGS) $(OF)e4sss \
		steadystate_solver_with_rev_string.d $(EILMER_FILES) \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(DLINKFLAGS)


eilmer-gpu-chem: main.d gpu_chem.d $(EILMER_FILES) \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) $(OPENCL_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) $(NM_FILES) $(UTIL_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' main.d > main_with_rev_string.d
	$(DMD) -version=gpu_chem $(DFLAGS) $(OF)eilmer-gpu-chem \
		main_with_rev_string.d $(EILMER_FILES) gpu_chem.d \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) $(OPENCL_FILES)\
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(DLINKFLAGS) -I$(OPENCL_DIR) -L-L/usr/lib64/nvidia -L-lOpenCL

prep-gpu-chem-kernel.py: utils/prep-gpu-chem-kernel.py
	cp utils/prep-gpu-chem-kernel.py .
	chmod +x prep-gpu-chem-kernel.py

alpha_qss_kernel.tmplt: utils/alpha_qss_kernel.tmplt
	cp utils/alpha_qss_kernel.tmplt .


