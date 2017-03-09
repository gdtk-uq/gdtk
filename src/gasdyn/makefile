# makefile for gasdyn
# Builds all of the demonstration programs.
#
# Peter J. September 2015
#          October 2016 bring into line with makefile from geom
#

PLATFORM ?= linux

DEMO_PROGRAMS := idealgasflow_demo luaidealgasflow_demo \
	gasflow_demo

include gasdyn_files.mk

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := $(LUA_DIR)/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

# For development work, we have put on -debug and -unittest.
DFLAGS := -w -debug -unittest -g -I..
OF := -of
DLINKFLAGS := -L-ldl

DFLAGS_LUA := -w -debug -unittest -g -I.. -I$(LUAD_DIR)
DLINKFLAGS_LUA := -L-L$(LIBLUAPATH) -L-llua $(DLINKFLAGS)

demo: $(DEMO_PROGRAMS)
	@echo "Demo codes built."

clean:
	- rm *.o
	- rm $(DEMO_PROGRAMS)
	- rm trace.log trace.def
	- rm thermo.inp thermo.lib thermo.out 
	- rm trans.inp trans.lib trans.out
	- rm tmp.inp tmp.out tmp.plt
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a

$(LIBLUA):
	cd $(LUA_DIR); make linux local

idealgasflow_demo: idealgasflow_demo.d idealgasflow.d $(NM_FILES) $(LIBLUA) \
	$(UTIL_FILES)
	dmd $(DFLAGS) $(OF)idealgasflow_demo \
		idealgasflow_demo.d idealgasflow.d $(NM_FILES) $(UTIL_FILES) \
		$(DLINKFLAGS_LUA)

luaidealgasflow_demo: luaidealgasflow_demo.d $(NM_FILES) $(LIBLUA) \
	$(UTIL_FILES) idealgasflow.d luaidealgasflow.d
	dmd $(DFLAGS_LUA) $(OF)luaidealgasflow_demo \
		luaidealgasflow_demo.d idealgasflow.d luaidealgasflow.d \
		$(NM_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA)

gasflow_demo: gasflow_demo.d gasflow.d $(NM_FILES) $(LIBLUA) \
	$(UTIL_FILES) $(GASDYN_FILES) $(GAS_FILES) $(KINETICS_FILES)
	dmd $(DFLAGS) $(OF)gasflow_demo \
		gasflow_demo.d gasflow.d $(NM_FILES) $(UTIL_FILES) \
		$(GAS_FILES) $(KINETICS_FILES) $(DLINKFLAGS_LUA)

#--------------------------------------------------------------------------
# Fortran bits -- we need a bit of a hack

$(GAS_DIR)/ideal_air_fortran.o: $(GAS_DIR)/ideal_air_fortran.f
	cd $(GAS_DIR); gfortran -c -ffree-form ideal_air_fortran.f
