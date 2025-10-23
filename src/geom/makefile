# makefile for the geometry package
# Builds all of the demonstration programs.

DMD ?= ldmd2

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR := ../../build
REPO_DIR := $(shell cd ../../; pwd)
ifneq ("$(wildcard $(REPO_DIR)/.git)","")
    # Repository exists
    REVISION_STRING := $(shell git rev-parse --short HEAD)
else
    REVISION_STRING := unknown
endif

DEMO_PROGRAMS := vector3_demo gpath_demo surface_demo volume_demo \
	univariatefunctions_demo \
	luageom_demo luagpath_demo luasurface_demo luavolume_demo \
	luaunifunction_demo \
	sgrid_demo luasgrid_demo usgrid_demo luausgrid_demo \
	paver_demo paver2d_demo svg_demo sketch_demo luasketch_demo \
	gpath_utils_demo

include geom_files.mk

NTYPES_DIR := ../ntypes
NTYPES_FILES := $(NTYPES_DIR)/complex.d

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib


DFLAGS := -w -debug -dip1008 -preview=in -g -I..
OF := -of
DLINKFLAGS :=
ifeq ($(findstring libplot,$(LIBPLOT_VERSION_STR)), libplot)
    DFLAGS := $(DFLAGS) -version=$(LIBPLOT_VERSION_STR)
    DLINKFLAGS := $(DLINKFLAGS) -L-lplot
endif
DLINKFLAGS := $(DLINKFLAGS) -L-ldl
ifeq ($(PLATFORM), macosx)
    DLINKFLAGS += -L-ld_classic
endif

DFLAGS_LUA := -w -dip1008 -g -I..
ifeq ($(findstring libplot,$(LIBPLOT_VERSION_STR)), libplot)
    DFLAGS_LUA := $(DFLAGS_LUA) -version=$(LIBPLOT_VERSION_STR)
endif

# ----------------------------------------------------------------------
# Here begins the list of targets, starting with the top-level actions.
#
# The default target is test-all, because that is most likely your interest
# if you are invoking this makefile from within the geom models directory.
# ----------------------------------------------------------------------

.PHONY: test-all test-real test-complex
test-all: test-real test-complex

test-real: test_runner_real
	@echo "Running real tests"
	- exec ./test_runner_real --module geom

test-complex: test_runner_complex
	@echo "Running complex tests"
	- exec ./test_runner_complex --module geom

# We need to pull the test_runner.d to the front to ensure no linker errors
test_runner_real: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES)
	@echo "Compiling test runner..."
	$(DMD) -of$@ -debug -g -w -dip1008 -unittest $(DLINKFLAGS) \
		../util/test_runner.d $(filter-out %/test_runner.d,$^)

test_runner_complex: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES)
	@echo "Compiling test runner with complex numbers..."
	$(DMD) -of$@ -debug -g -w -dip1008 -unittest -version=complex_numbers $(DLINKFLAGS) \
		../util/test_runner.d $(filter-out %/test_runner.d,$^)

play:
	@echo "Found $(LIBPLOT)"

demo: $(DEMO_PROGRAMS)
	@echo "Demo codes built."

install: foamMesh grid2pato $(GEOM_LUA_FILES)
	- mkdir -p $(INSTALL_DIR)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share
	cp foamMesh $(BUILD_DIR)/bin/
	cp $(OPENFOAM_DIR)/foam-mesh.lua $(BUILD_DIR)/bin/
	cp $(OPENFOAM_DIR)/openfoam.lua $(BUILD_DIR)/lib/
	cp -r $(OPENFOAM_DIR)/foamMesh-templates $(BUILD_DIR)/share/
	cp $(OPENFOAM_DIR)/grid2pato.lua $(BUILD_DIR)/bin/
	cp grid2pato $(BUILD_DIR)/bin/
	cp -r $(OPENFOAM_DIR)/pato-templates $(BUILD_DIR)/share/
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

clean:
	- rm -f *.o elements/*.o gpath/*.o grid/*.o luawrap/*.o misc/*.o surface/*.o volume/*.o
	- rm -f *.obj elements/*.obj gpath/*.obj grid/*.obj luawrap/*.obj misc/*.obj surface/*.obj volume/*.obj
	- rm -f $(DEMO_PROGRAMS)
	- rm -f test_runner_real test_runner_complex
	- rm -f gpgrid*.vtk test*.vtk test*.stl test*.gz test*.svg test*.su2 test*.bin
	- rm -rf test_openFoam
	- rm $(BUILD_DIR)/lib/openfoam.lua
	- cd $(LUA_DIR); make clean

all-clean:
	- make clean
	- rm -rf foamMesh $(OPENFOAM_DIR)/foamMesh_with_rev_string.d
	- rm -rf grid2pato $(OPENFOAM_DIR)/grid2pato_with_rev_string.d
	- rm -rf $(BUILD_DIR)

# Specific targets, for individual artefacts.

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

foamMesh: $(OPENFOAM_DIR)/foamMesh.d $(LIBLUA) $(GEOM_FILES) $(GEOM_LUA_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' $(OPENFOAM_DIR)/foamMesh.d > $(OPENFOAM_DIR)/foamMesh_with_rev_string.d

	$(DMD) $(DFLAGS_LUA) $(OPENFOAM_DIR)/foamMesh_with_rev_string.d $(OF)foamMesh \
	$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

grid2pato: $(OPENFOAM_DIR)/grid2pato.d $(LIBLUA) $(GEOM_FILES) $(GEOM_LUA_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' $(OPENFOAM_DIR)/grid2pato.d > $(OPENFOAM_DIR)/grid2pato_with_rev_string.d

	$(DMD) $(DFLAGS_LUA) $(OPENFOAM_DIR)/grid2pato_with_rev_string.d $(OF)grid2pato \
	$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)


# The original demo artefacts.

sketch_demo: misc/sketch_demo.d $(GEOM_D_FILES) $(NM_FILES) $(NTYPES_FILES)
	$(DMD) $(DFLAGS) misc/sketch_demo.d $(OF)sketch_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(DLINKFLAGS)

svg_demo: misc/svg_demo.d misc/svg.d
	$(DMD) $(DFLAGS) misc/svg_demo.d $(OF)svg_demo misc/svg.d $(DLINKFLAGS)

vector3_demo: elements/vector3_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES)
	$(DMD) $(DFLAGS) elements/vector3_demo.d $(OF)vector3_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(DLINKFLAGS)

gpath_demo: gpath/gpath_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES)
	$(DMD) $(DFLAGS) gpath/gpath_demo.d $(OF)gpath_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(DLINKFLAGS)

gpath_utils_demo: gpath/gpath_utils_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES)
	$(DMD) $(DFLAGS) gpath/gpath_utils_demo.d $(OF)gpath_utils_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(DLINKFLAGS)

surface_demo: surface/surface_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES)
	$(DMD) $(DFLAGS) surface/surface_demo.d $(OF)surface_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(DLINKFLAGS)

volume_demo: volume/volume_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES)
	$(DMD) $(DFLAGS) volume/volume_demo.d $(OF)volume_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(DLINKFLAGS)

univariatefunctions_demo: misc/univariatefunctions_demo.d misc/univariatefunctions.d
	$(DMD) $(DFLAGS) misc/univariatefunctions_demo.d $(OF)univariatefunctions_demo \
		misc/univariatefunctions.d $(DLINKFLAGS)

sgrid_demo: grid/sgrid_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES)
	$(DMD) $(DFLAGS) grid/sgrid_demo.d $(OF)sgrid_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(DLINKFLAGS)

usgrid_demo: grid/usgrid_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES)
	$(DMD) $(DFLAGS) grid/usgrid_demo.d $(OF)usgrid_demo \
		 $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(DLINKFLAGS)

luageom_demo: luawrap/luageom_demo.d $(GEOM_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA)
	$(DMD) $(DFLAGS_LUA) luawrap/luageom_demo.d $(OF)luageom_demo \
		$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

luagpath_demo: luawrap/luagpath_demo.d $(GEOM_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA)
	$(DMD) $(DFLAGS_LUA) luawrap/luagpath_demo.d $(OF)luagpath_demo \
		$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

luasurface_demo: luawrap/luasurface_demo.d $(GEOM_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA)
	$(DMD) $(DFLAGS_LUA) luawrap/luasurface_demo.d $(OF)luasurface_demo \
		$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

luasketch_demo: luawrap/luasketch_demo.d $(GEOM_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA)
	$(DMD) $(DFLAGS_LUA) luawrap/luasketch_demo.d $(OF)luasketch_demo \
		$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

luavolume_demo: luawrap/luavolume_demo.d $(GEOM_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA)
	$(DMD) $(DFLAGS_LUA) luawrap/luavolume_demo.d $(OF)luavolume_demo \
		$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

luaunifunction_demo: luawrap/luaunifunction_demo.d luawrap/luaunifunction.d \
	misc/univariatefunctions.d $(UTIL_FILES) $(LIBLUA)
	$(DMD) $(DFLAGS_LUA) luawrap/luaunifunction_demo.d $(OF)luaunifunction_demo \
		luawrap/luaunifunction.d misc/univariatefunctions.d \
		$(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

luasgrid_demo: luawrap/luasgrid_demo.d $(GEOM_FILES) $(NM_FILES) $(NTYPES_FILES) $(GZIP_FILES) $(UTIL_FILES)
	$(DMD) $(DFLAGS_LUA) luawrap/luasgrid_demo.d $(OF)luasgrid_demo \
		$(GEOM_FILES) $(NM_FILES) $(NTYPES_FILES) $(GZIP_FILES) $(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

luausgrid_demo: luawrap/luausgrid_demo.d $(GEOM_FILES) $(NM_FILES) $(NTYPES_FILES) $(GZIP_FILES) $(UTIL_FILES)
	$(DMD) $(DFLAGS_LUA) luawrap/luausgrid_demo.d $(OF)luausgrid_demo \
		$(GEOM_FILES) $(NM_FILES) $(NTYPES_FILES) $(GZIP_FILES) $(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

paver_demo: grid/paver_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES)
	$(DMD) $(DFLAGS) grid/paver_demo.d $(OF)paver_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(DLINKFLAGS)

paver2d_demo: grid/paver2d_demo.d $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES)
	$(DMD) $(DFLAGS) grid/paver2d_demo.d $(OF)paver2d_demo \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(DLINKFLAGS)
