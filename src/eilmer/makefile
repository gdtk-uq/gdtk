# makefile for Eilmer4
# Builds main simulation program by default.

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd

PROGRAMS := e4shared

DEMO_PROGRAMS := globalconfig_demo flowstate_demo fv_demo block_demo \
	luaflowstate_demo luaglobalconfig_demo

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

GAS_DIR := ../gas

GAS_MODEL_FILES := \
	$(GAS_DIR)/package.d \
	$(GAS_DIR)/gas_model.d \
	$(GAS_DIR)/gas_model_util.d \
	$(GAS_DIR)/ideal_gas.d \
	$(GAS_DIR)/physical_constants.d \
	$(GAS_DIR)/therm_perf_gas.d \
        $(GAS_DIR)/very_viscous_air.d\
	$(GAS_DIR)/CO2Gas.d

THERMO_FILES := \
	$(GAS_DIR)/thermo/cea_thermo_curves.d \
	$(GAS_DIR)/thermo/evt_eos.d \
	$(GAS_DIR)/thermo/perf_gas_mix_eos.d \
	$(GAS_DIR)/thermo/pvt_eos.d \
	$(GAS_DIR)/thermo/therm_perf_gas_mix_eos.d

DIFFUSION_FILES := \
	$(GAS_DIR)/diffusion/cea_therm_cond.d \
	$(GAS_DIR)/diffusion/cea_viscosity.d \
	$(GAS_DIR)/diffusion/sutherland_therm_cond.d \
	$(GAS_DIR)/diffusion/sutherland_viscosity.d \
	$(GAS_DIR)/diffusion/therm_cond.d \
	$(GAS_DIR)/diffusion/viscosity.d \
	$(GAS_DIR)/diffusion/wilke_mixing_therm_cond.d \
	$(GAS_DIR)/diffusion/wilke_mixing_viscosity.d

GAS_FILES := $(GAS_MODEL_FILES) $(THERMO_FILES) $(DIFFUSION_FILES)

KINETICS_DIR := ../kinetics

KINETICS_FILES := \
	$(KINETICS_DIR)/package.d \
	$(KINETICS_DIR)/chemistry_update.d \
	$(KINETICS_DIR)/rate_constant.d \
	$(KINETICS_DIR)/reaction.d \
	$(KINETICS_DIR)/reaction_mechanism.d

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua.d \
	$(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

NM_DIR := ../nm
NM_FILES := $(NM_DIR)/ridder.d \
	$(NM_DIR)/bbla.d

GEOM_DIR := ../geom
GEOM_FILES := $(GEOM_DIR)/geom.d \
	$(GEOM_DIR)/univariatefunctions.d \
	$(GEOM_DIR)/luaunifunction.d \
	$(GEOM_DIR)/gpath.d \
	$(GEOM_DIR)/surface.d \
	$(GEOM_DIR)/volume.d \
	$(GEOM_DIR)/luageom.d \
	$(GEOM_DIR)/luagpath.d \
	$(GEOM_DIR)/luasurface.d \
	$(GEOM_DIR)/luavolume.d \
	$(GEOM_DIR)/sgrid.d \
	$(GEOM_DIR)/luasgrid.d

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := $(LUA_DIR)/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib

# The install destination the same as Eilmer3, for now.
INSTALL_DIR ?= $(HOME)/e3bin

REVISION_STRING := $(shell hg identify --id --num --branch --tags)

ifeq ($(DMD), dmd)
    # (1) For debug build.  ([TODO] fix -unittest build.)
    #     DFLAGS := -debug -gc -w
    # (2) For profiling the debugged code.
    #     DFLAGS := -profile -vgc -w -O -release -inline -boundscheck=off
    # (3) For releasing a fast code.
    #     DFLAGS := -w -O -release -inline -boundscheck=off
    DFLAGS := -w -O -release -inline -boundscheck=off \
	-I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GZIP_DIR)
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -O -release -inline -boundscheck=off \
	-I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GZIP_DIR)
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), gdc)
    # For debug: DFLAGS := -Wall -Og -g
    DFLAGS := -Wall -O2 -frelease -fno-debug \
	-I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GZIP_DIR)
    OF := -o
    DLINKFLAGS := -L$(LIBLUAPATH) $(LIBLUA) -ldl
endif

default: $(PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Eilmer4 simulation code built."

install: $(PROGRAMS) prep.lua prep-gas prep-chem
	@echo "Installing to $(INSTALL_DIR)"
	cp $(PROGRAMS) $(INSTALL_DIR)
	cp prep.lua $(INSTALL_DIR)

demo: $(DEMO_PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Demo codes built."

clean:
	- rm *.o
	- rm $(DEMO_PROGRAMS) $(PROGRAMS)
	- rm test-grid.txt.gz test-flow.txt.gz trace.log trace.def
	- rm -r plot/ flow/ grid/
	- rm main_with_rev_string.d
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a

prep-gas:
	cd $(GAS_DIR); make INSTALL_DIR=$(INSTALL_DIR) install-prep-gas 

prep-chem:
	cd $(KINETICS_DIR); make INSTALL_DIR=$(INSTALL_DIR) install-prep-chem

$(LIBLUA):
	cd $(LUA_DIR); make linux local

globalconfig_demo: globalconfig_demo.d globalconfig.d fvcore.d \
	$(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) \
	$(KINETICS_FILES) $(GZIP_FILES)
	$(DMD) $(DFLAGS) $(OF)globalconfig_demo \
		globalconfig_demo.d globalconfig.d fvcore.d \
		$(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)  

luaglobalconfig_demo: luaglobalconfig_demo.d luaglobalconfig.d globalconfig.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d conservedquantities.d \
	flowstate.d json_helper.d \
	$(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) \
	$(KINETICS_FILES) $(GZIP_FILES)
	$(DMD) $(DFLAGS) $(OF)luaglobalconfig_demo \
		luaglobalconfig_demo.d luaglobalconfig.d globalconfig.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d conservedquantities.d \
		flowstate.d json_helper.d \
		$(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(DLINKFLAGS)  

flowstate_demo: flowstate_demo.d flowstate.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d conservedquantities.d \
	globalconfig.d json_helper.d \
	$(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) \
	$(KINETICS_FILES) $(GZIP_FILES)
	$(DMD) $(DFLAGS) $(OF)flowstate_demo flowstate_demo.d \
		fvcore.d flowstate.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d globalconfig.d json_helper.d \
		$(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(DLINKFLAGS)

luaflowstate_demo: luaflowstate_demo.d luaflowstate.d flowstate.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d conservedquantities.d \
	globalconfig.d luaglobalconfig.d json_helper.d \
	$(GEOM_FILES) $(BC_FILES) $(SOLID_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) \
	$(KINETICS_FILES) $(GZIP_FILES)
	$(DMD) $(DFLAGS) $(OF)luaflowstate_demo \
		luaflowstate_demo.d luaflowstate.d flowstate.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d conservedquantities.d \
		globalconfig.d luaglobalconfig.d json_helper.d \
		$(GEOM_FILES) $(BC_FILES) $(SOLID_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(DLINKFLAGS)

fv_demo: fv_demo.d fvcore.d fvcell.d fvvertex.d fvinterface.d onedinterp.d fluxcalc.d \
	conservedquantities.d flowstate.d luaflowstate.d \
	globalconfig.d globaldata.d \
	block.d sblock.d ublock.d \
	json_helper.d lua_helper.d user_defined_source_terms.d \
	simcore.d fileutil.d readconfig.d luaglobalconfig.d \
	$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) 
	$(DMD) $(DFLAGS) $(OF)fv_demo \
		fv_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d flowstate.d luaflowstate.d onedinterp.d fluxcalc.d \
		globalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		json_helper.d lua_helper.d user_defined_source_terms.d \
		simcore.d fileutil.d readconfig.d luaglobalconfig.d \
		$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(DLINKFLAGS)

block_demo: block_demo.d \
	block.d sblock.d ublock.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d \
	fluxcalc.d onedinterp.d viscousflux.d \
	conservedquantities.d flowstate.d luaflowstate.d \
	globalconfig.d globaldata.d \
	block.d sblock.d ublock.d \
	json_helper.d lua_helper.d user_defined_source_terms.d \
	simcore.d fileutil.d readconfig.d luaglobalconfig.d \
	$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) 
	$(DMD) $(DFLAGS) $(OF)block_demo \
		block_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		fluxcalc.d onedinterp.d viscousflux.d \
		conservedquantities.d flowstate.d luaflowstate.d \
		globalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		json_helper.d lua_helper.d user_defined_source_terms.d \
		simcore.d fileutil.d readconfig.d luaglobalconfig.d \
		$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(DLINKFLAGS)

e4shared: main.d simcore.d fileutil.d readconfig.d json_helper.d \
	block.d sblock.d ublock.d \
	fluxcalc.d onedinterp.d viscousflux.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d \
	conservedquantities.d flowstate.d luaflowstate.d \
	globalconfig.d luaglobalconfig.d globaldata.d \
	block.d sblock.d ublock.d \
	lua_helper.d user_defined_source_terms.d gas_solid_interface.d \
	postprocess.d flowsolution.d \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) 
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' main.d > main_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4shared \
		main_with_rev_string.d simcore.d fileutil.d readconfig.d json_helper.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		fluxcalc.d onedinterp.d viscousflux.d \
		conservedquantities.d flowstate.d luaflowstate.d \
		globalconfig.d luaglobalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		lua_helper.d user_defined_source_terms.d gas_solid_interface.d \
		postprocess.d flowsolution.d \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(DLINKFLAGS)

compute-error: compute_error.d simcore.d fileutil.d readconfig.d json_helper.d \
	block.d sblock.d ublock.d \
	fluxcalc.d onedinterp.d viscousflux.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d \
	conservedquantities.d flowstate.d luaflowstate.d \
	globalconfig.d luaglobalconfig.d globaldata.d \
	block.d sblock.d ublock.d \
	lua_helper.d user_defined_source_terms.d gas_solid_interface.d \
	$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(UTIL_FILES) $(NM_FILES) 
	$(DMD) $(DFLAGS) $(OF)compute-error \
		compute_error.d simcore.d fileutil.d readconfig.d json_helper.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		fluxcalc.d onedinterp.d viscousflux.d \
		conservedquantities.d flowstate.d luaflowstate.d \
		globalconfig.d luaglobalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		lua_helper.d user_defined_source_terms.d gas_solid_interface.d \
		$(BC_FILES) $(SOLID_FILES) $(GEOM_FILES) $(GAS_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(KINETICS_FILES) $(DLINKFLAGS)
