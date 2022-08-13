# makefile for gasdyn
# Builds all of the demonstration programs.
#
# Peter J. September 2015
#          October 2016 bring into line with makefile from geom
#

PLATFORM ?= linux

DEMO_PROGRAMS := idealgasflow_demo luaidealgasflow_demo \
	gasflow_demo luagasflow_demo

include gasdyn_files.mk

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib

CEQ_DIR := ../extern/ceq/source
LIBCEQ := $(CEQ_DIR)/libceq.a
include $(CEQ_DIR)/ceq_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk
LIBGASF := $(GAS_DIR)/libgasf.a

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

# For development work, we have put on -debug.
DFLAGS := -w -debug -g -dip1008 -I..
OF := -of
DLINKFLAGS := -L-ldl -L-lgfortran

DFLAGS_LUA := -w -debug -g -dip1008 -I.. -I$(LUAD_DIR)
DLINKFLAGS_LUA := -L-L$(LIBLUAPATH) $(DLINKFLAGS)

demo: $(DEMO_PROGRAMS)
	@echo "Demo codes built."

clean:
	- rm *.o *.obj
	- rm $(DEMO_PROGRAMS)
	- rm trace.log trace.def
	- rm thermo.inp thermo.lib thermo.out
	- rm trans.inp trans.lib trans.out
	- rm tmp.inp tmp.out tmp.plt
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(CEQ_DIR); make clean

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local

$(LIBGASF):
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) libgasf.a

$(LIBCEQ):
	cd $(CEQ_DIR); make

idealgasflow_demo: idealgasflow_demo.d idealgasflow.d $(NM_FILES) $(LIBLUA) \
	$(UTIL_FILES)
	dmd $(DFLAGS) $(OF)idealgasflow_demo \
		idealgasflow_demo.d idealgasflow.d gasflowexception.d \
		$(NM_FILES) $(LIBLUA) $(UTIL_FILES) $(DLINKFLAGS_LUA)

luaidealgasflow_demo: luaidealgasflow_demo.d $(NM_FILES) $(LIBLUA) \
	$(UTIL_FILES) idealgasflow.d luaidealgasflow.d
	dmd $(DFLAGS_LUA) $(OF)luaidealgasflow_demo \
		luaidealgasflow_demo.d idealgasflow.d luaidealgasflow.d \
		gasflowexception.d \
		$(NM_FILES) $(LIBLUA) $(UTIL_FILES) $(DLINKFLAGS_LUA)

gasflow_demo: gasflow_demo.d gasflow.d idealgasflow.d $(NM_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) \
	$(UTIL_FILES) $(GASDYN_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES)
	dmd $(DFLAGS) $(OF)gasflow_demo \
		gasflow_demo.d gasflow.d idealgasflow.d gasflowexception.d \
		$(NM_FILES) $(UTIL_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(KINETICS_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(DLINKFLAGS_LUA)

luagasflow_demo: luagasflow_demo.d luagasflow.d luaidealgasflow.d \
	gasflow.d idealgasflow.d \
	$(NM_FILES) $(LIBLUA) $(LIBCEQ) $(LIBGASF) \
	$(UTIL_FILES) $(GASDYN_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(GAS_LUA_FILES) \
	$(KINETICS_FILES) $(KINETICS_LUA_FILES)
	dmd $(DFLAGS) $(OF)luagasflow_demo \
		luagasflow_demo.d luagasflow.d luaidealgasflow.d \
		gasflow.d idealgasflow.d gasflowexception.d \$(NM_FILES) $(UTIL_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GAS_LUA_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) $(DLINKFLAGS_LUA)
