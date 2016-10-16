# makefile for gasdyn
# Builds all of the demonstration programs.
#
# Peter J. September 2015
#          October 2016 bring into line with makefile from geom
#

PLATFORM ?= linux

DEMO_PROGRAMS := idealgasflow_demo luaidealgasflow_demo

include gasdyn_files.mk

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua.d \
	$(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

NM_DIR := ../nm
NM_FILES := $(NM_DIR)/ridder.d $(NM_DIR)/bracketing.d \
	$(NM_DIR)/linesearch.d $(NM_DIR)/bbla.d

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := $(LUA_DIR)/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := $(wildcard $(LUAD_DIR)/luad/*.d) \
	$(wildcard $(LUAD_DIR)/luad/c/*.d) \
	$(wildcard $(LUAD_DIR)/luad/conversions/*.d)

DFLAGS := -w -unittest -g -I..
OF := -of
DLINKFLAGS := -L-ldl

DFLAGS_LUA := -w -g -I.. -I$(LUAD_DIR)
DLINKFLAGS_LUA := -L-L$(LIBLUAPATH) -L-llua $(DLINKFLAGS)

demo: $(DEMO_PROGRAMS)
	@echo "Demo codes built."

clean:
	- rm *.o
	- rm $(DEMO_PROGRAMS)
	- rm trace.log trace.def
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a

$(LIBLUA):
	cd $(LUA_DIR); make linux local

idealgasflow_demo: idealgasflow_demo.d idealgasflow.d $(NM_FILES) $(LIBLUA) \
	$(UTIL_FILES) $(GASDYN_FILES)
	dmd $(DFLAGS) $(OF)idealgasflow_demo \
		idealgasflow_demo.d idealgasflow.d $(NM_FILES) $(DLINKFLAGS)

luaidealgasflow_demo: luaidealgasflow_demo.d $(NM_FILES) $(LIBLUA) \
	$(UTIL_FILES) $(GASDYN_FILES)
	dmd $(DFLAGS_LUA) $(OF)luaidealgasflow_demo \
		luaidealgasflow_demo.d $(GASDYN_FILES) $(NM_FILES) \
		$(UTIL_FILES) $(DLINKFLAGS_LUA)
