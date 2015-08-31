# makefile for the geometry package
# Builds all of the demonstration programs.

DEMO_PROGRAMS := geom_demo gpath_demo surface_demo volume_demo \
	univariatefunctions_demo \
	luageom_demo luagpath_demo luasurface_demo luavolume_demo \
	luaunifunction_demo \
	sgrid_demo luasgrid_demo

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua.d \
	$(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

NM_DIR := ../nm
NM_FILES := $(NM_DIR)/ridder.d \
	$(NM_DIR)/bbla.d

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
	echo "Demo codes built."

clean:
	- rm *.o
	- rm $(DEMO_PROGRAMS)
	- cd $(LUA_DIR); make clean

$(LIBLUA):
	cd $(LUA_DIR); make linux local

geom_demo: geom_demo.d geom.d
	dmd $(DFLAGS) geom_demo.d $(OF)geom_demo \
		geom.d $(DLINKFLAGS)

gpath_demo: gpath_demo.d gpath.d geom.d $(NM_FILES)
	dmd $(DFLAGS) gpath_demo.d $(OF)gpath_demo \
		gpath.d geom.d $(NM_FILES) $(DLINKFLAGS)

surface_demo: surface_demo.d surface.d gpath.d geom.d $(NM_FILES)
	dmd $(DFLAGS) surface_demo.d $(OF)surface_demo \
		surface.d gpath.d geom.d $(NM_FILES) $(DLINKFLAGS)

volume_demo: volume_demo.d volume.d surface.d gpath.d geom.d $(NM_FILES)
	dmd $(DFLAGS) volume_demo.d $(OF)volume_demo \
		volume.d surface.d gpath.d geom.d $(NM_FILES) $(DLINKFLAGS)

univariatefunctions_demo: univariatefunctions_demo.d univariatefunctions.d
	dmd $(DFLAGS) univariatefunctions_demo.d $(OF)univariatefunctions_demo \
		univariatefunctions.d $(DLINKFLAGS)

sgrid_demo: sgrid_demo.d sgrid.d volume.d surface.d gpath.d geom.d \
	univariatefunctions.d \
	$(GZIP_FILES) $(NM_FILES)
	dmd $(DFLAGS) sgrid_demo.d $(OF)sgrid_demo \
		sgrid.d volume.d surface.d gpath.d geom.d \
		univariatefunctions.d \
		$(GZIP_FILES) $(NM_FILES) $(DLINKFLAGS)

luageom_demo: luageom_demo.d luageom.d geom.d $(LIBLUA)
	dmd $(DFLAGS_LUA) luageom_demo.d $(OF)luageom_demo \
		luageom.d geom.d $(UTIL_FILES) $(DLINKFLAGS_LUA)

luagpath_demo: luagpath_demo.d luagpath.d luageom.d geom.d gpath.d \
	$(NM_FILES) $(UTIL_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luagpath_demo.d $(OF)luagpath_demo \
		luagpath.d luageom.d geom.d gpath.d \
		$(NM_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA)

luasurface_demo: luasurface_demo.d luasurface.d luagpath.d luageom.d \
	geom.d gpath.d surface.d \
	$(NM_FILES) $(UTIL_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luasurface_demo.d $(OF)luasurface_demo \
		luasurface.d luagpath.d luageom.d \
		surface.d geom.d gpath.d \
		$(NM_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA)

luavolume_demo: luavolume_demo.d luavolume.d luasurface.d luagpath.d luageom.d \
	geom.d gpath.d surface.d volume.d \
	$(NM_FILES) $(UTIL_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luavolume_demo.d $(OF)luavolume_demo \
		luavolume.d luasurface.d luagpath.d luageom.d \
		volume.d surface.d geom.d gpath.d \
		$(NM_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA)

luaunifunction_demo: luaunifunction_demo.d luaunifunction.d univariatefunctions.d \
	$(UTIL_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luaunifunction_demo.d $(OF)luaunifunction_demo \
		luaunifunction.d univariatefunctions.d \
		$(UTIL_FILES) $(DLINKFLAGS_LUA)

luasgrid_demo: luasgrid_demo.d sgrid.d volume.d surface.d gpath.d geom.d \
	univariatefunctions.d \
	luasgrid.d luavolume.d luasurface.d luagpath.d luageom.d \
	luaunifunction.d \
	$(NM_FILES) $(GZIP_FILES) $(UTIL_FILES) 
	dmd $(DFLAGS_LUA) luasgrid_demo.d $(OF)luasgrid_demo \
		sgrid.d volume.d surface.d gpath.d geom.d \
		univariatefunctions.d \
		luasgrid.d luavolume.d luasurface.d luagpath.d luageom.d \
		luaunifunction.d \
		$(NM_FILES) $(GZIP_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA)

