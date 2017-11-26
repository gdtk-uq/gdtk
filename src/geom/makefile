# makefile for the geometry package
# Builds all of the demonstration programs.

DMD ?= dmd
PLATFORM ?= linux
INSTALL_DIR ?= $(HOME)/dgdinst
BUILD_DIR := ../../build
REVISION_STRING := $(shell hg identify --id --num --branch --tags)

DEMO_PROGRAMS := vector3_demo gpath_demo surface_demo volume_demo \
	univariatefunctions_demo \
	luageom_demo luagpath_demo luasurface_demo luavolume_demo \
	luaunifunction_demo \
	sgrid_demo luasgrid_demo usgrid_demo luausgrid_demo \
	paver_demo svg_demo sketch_demo luasketch_demo

include geom_files.mk

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
	- rm gpgrid*.vtk test*.vtk test*.gz test*.svg test*.su2 test*.bin
	- rm -r test_openFoam
	- cd $(LUA_DIR); make clean

install: foamMesh foam-mesh.lua 
	- mkdir -p $(INSTALL_DIR)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share
	cp foamMesh $(BUILD_DIR)/bin/
	cp foam-mesh.lua $(BUILD_DIR)/bin/
	cp -r foamMesh-templates $(BUILD_DIR)/share/
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)


$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local

foamMesh: foamMesh.d $(LIBLUA) $(GEOM_FILES) $(GEOM_LUA_FILES) $(GZIP_FILES) $(NM_FILES) $(UTIL_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' foamMesh.d > foamMesh_with_rev_string.d

	dmd $(DFLAGS_LUA) foamMesh_with_rev_string.d $(OF)foamMesh \
	$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA) 

sketch_demo: sketch_demo.d $(GEOM_D_FILES) $(NM_FILES)
	dmd $(DFLAGS) sketch_demo.d $(OF)sketch_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(DLINKFLAGS)

svg_demo: svg_demo.d svg.d
	dmd $(DFLAGS) svg_demo.d $(OF)svg_demo svg.d $(DLINKFLAGS)

vector3_demo: elements/vector3_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES)
	dmd $(DFLAGS) elements/vector3_demo.d $(OF)vector3_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(DLINKFLAGS)

gpath_demo: gpath/gpath_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES)
	dmd $(DFLAGS) gpath/gpath_demo.d $(OF)gpath_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(DLINKFLAGS)

surface_demo: surface_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES)
	dmd $(DFLAGS) surface_demo.d $(OF)surface_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(DLINKFLAGS)

volume_demo: volume_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES)
	dmd $(DFLAGS) volume_demo.d $(OF)volume_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(DLINKFLAGS)

univariatefunctions_demo: univariatefunctions_demo.d univariatefunctions.d
	dmd $(DFLAGS) univariatefunctions_demo.d $(OF)univariatefunctions_demo \
		univariatefunctions.d $(DLINKFLAGS)

sgrid_demo: sgrid_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES)
	dmd $(DFLAGS) sgrid_demo.d $(OF)sgrid_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(DLINKFLAGS)

usgrid_demo: usgrid_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES)
	dmd $(DFLAGS) usgrid_demo.d $(OF)usgrid_demo \
		 $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(DLINKFLAGS)

luageom_demo: luageom_demo.d $(GEOM_FILES) $(NM_FILES) $(UTIL_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luageom_demo.d $(OF)luageom_demo \
		$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) \
		$(UTIL_FILES) $(DLINKFLAGS_LUA)

luagpath_demo: luagpath_demo.d $(GEOM_FILES) $(NM_FILES) $(UTIL_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luagpath_demo.d $(OF)luagpath_demo \
		$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) \
		$(UTIL_FILES) $(DLINKFLAGS_LUA)

luasurface_demo: luasurface_demo.d $(GEOM_FILES) $(NM_FILES) $(UTIL_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luasurface_demo.d $(OF)luasurface_demo \
		$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA)

luasketch_demo: luasketch_demo.d $(GEOM_FILES) $(NM_FILES) $(UTIL_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luasketch_demo.d $(OF)luasketch_demo \
		$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA)

luavolume_demo: luavolume_demo.d $(GEOM_FILES) $(NM_FILES) $(UTIL_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luavolume_demo.d $(OF)luavolume_demo \
		$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA)

luaunifunction_demo: luaunifunction_demo.d luaunifunction.d univariatefunctions.d \
	$(UTIL_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luaunifunction_demo.d $(OF)luaunifunction_demo \
		luaunifunction.d univariatefunctions.d \
		$(UTIL_FILES) $(DLINKFLAGS_LUA)

luasgrid_demo: luasgrid_demo.d $(GEOM_FILES) $(NM_FILES) $(GZIP_FILES) $(UTIL_FILES) 
	dmd $(DFLAGS_LUA) luasgrid_demo.d $(OF)luasgrid_demo \
		$(GEOM_FILES) $(NM_FILES) $(GZIP_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA)

luausgrid_demo: luausgrid_demo.d $(GEOM_FILES) $(NM_FILES) $(GZIP_FILES) $(UTIL_FILES) 
	dmd $(DFLAGS_LUA) luausgrid_demo.d $(OF)luausgrid_demo \
		$(GEOM_FILES) $(NM_FILES) $(GZIP_FILES) $(UTIL_FILES) $(DLINKFLAGS_LUA)

paver_demo: paver_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES)
	dmd $(DFLAGS) paver_demo.d $(OF)paver_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(DLINKFLAGS)

