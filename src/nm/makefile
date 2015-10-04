# makefile for the Numerical Methods package
# Builds all of the demonstration programs
# includes make for build_tree which demonstrates the tree_patch

DEMO_PROGRAMS := ridder_demo linesearch_demo nelmin_demo rungekutta_demo \
	newtoncotes_demo gaussquad_demo bbla_demo build_tree univariate_lut_demo

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua.d \
	$(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

NM_DIR := .
NM_FILES := $(NM_DIR)/ridder.d\
			$(NM_DIR)/tree_patch.d\
			$(NM_DIR)/univariate_lut.d 


LUA := ../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib



GAS_DIR := ../gas
GAS_MODEL_FILES := \
	$(GAS_DIR)/package.d \
	$(GAS_DIR)/gas_model.d \
	$(GAS_DIR)/gas_model_util.d \
	$(GAS_DIR)/ideal_gas.d \
	$(GAS_DIR)/physical_constants.d \
	$(GAS_DIR)/therm_perf_gas.d \
	$(GAS_DIR)/very_viscous_air.d\
	$(GAS_DIR)/CO2Gas.d\
	$(GAS_DIR)/CO2GasSW.d\
	$(GAS_DIR)/SF6Virial.d

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

DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl

libgas.a : $(GAS_FILES) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -lib $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(OF)$(LIBNAME) $(DFLAGS)

$(LIBLUA):
	cd $(LUA); make linux local

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

ridder_demo: ridder_demo.d ridder.d
	dmd -unittest -g ridder_demo.d

# Don't know why we have to include linesearch.d on the command line.
linesearch_demo: linesearch_demo.d linesearch.d
	dmd -unittest -g linesearch_demo.d linesearch.d

nelmin_demo: nelmin_demo.d nelmin.d
	dmd -unittest -g nelmin_demo.d nelmin.d

rungekutta_demo: rungekutta_demo.d rungekutta.d
	dmd -unittest -g rungekutta_demo.d rungekutta.d

newtoncotes_demo: newtoncotes_demo.d newtoncotes.d
	dmd -unittest -g newtoncotes_demo.d newtoncotes.d

gaussquad_demo: gaussquad_demo.d gaussquad.d
	dmd -unittest -g gaussquad_demo.d gaussquad.d

bbla_demo: bbla_demo.d bbla.d
	dmd -unittest -g bbla_demo.d bbla.d

build_tree: build_tree.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	dmd build_tree.d $(GAS_FILES) $(UTIL_FILES) $(NM_FILES)\
		$(DLINKFLAGS)

univariate_lut_demo:  univariate_lut_demo.d $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	dmd -debug univariate_lut_demo.d $(GAS_FILES) $(UTIL_FILES) $(NM_FILES)\
		$(DLINKFLAGS)
