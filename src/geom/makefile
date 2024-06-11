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

TEST_PROGRAMS := vector3_test vector3_complex_test \
	projection_test projection_complex_test \
	properties_test properties_complex_test \
	line_test line_complex_test \
	arc_test arc_complex_test \
	bezier_test bezier_complex_test \
	helix_test helix_complex_test \
	nurbs_utils_test nurbs_test \
	polyline_test polyline_complex_test \
	xspline_test \
	xsplinelsq_test \
	svgpath_test \
	coonspatch_test coonspatch_complex_test \
	aopatch_test aopatch_complex_test \
	gmopatch_test gmopatch_complex_test \
	bezierpatch_test bezierpatch_complex_test \
	controlpointpatch_test \
	nurbssurface_test \
	tfivolume_test tfivolume_complex_test \
	sweptsurfacevolume_test sweptsurfacevolume_complex_test \
	wedgevolume_test wedgevolume_complex_test \
	slabvolume_test slabvolume_complex_test \
	nurbsvolume_test \
	twosurfacevolume_test twosurfacevolume_complex_test \
	sgrid_test sgrid_complex_test \
	usgrid_test usgrid_complex_test \
	univariatefunctions_test svg_test \
	beziertrianglepatch_test \
	spherepatch_test spherepatch_complex_test \
	cubepatch_test cubepatch_complex_test

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
UTIL_FILES := $(UTIL_DIR)/lua.d \
	$(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

NM_DIR := ../nm
NM_FILES := $(NM_DIR)/ridder.d \
	$(NM_DIR)/bbla.d \
	$(NM_DIR)/linesearch.d \
	$(NM_DIR)/spline.d \
	$(NM_DIR)/splinelsq.d \
	$(NM_DIR)/nelmin.d \
	$(NM_DIR)/number.d

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
# The default target is test, because that is most likely your interest
# if you are invoking this makefile from within the gas models directory.
# ----------------------------------------------------------------------
play:
	@echo "Found $(LIBPLOT)"

test: $(TEST_PROGRAMS)
	tclsh geom-package-test.tcl

demo: $(DEMO_PROGRAMS)
	@echo "Demo codes built."

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

clean:
	- rm -f *.o elements/*.o gpath/*.o grid/*.o luawrap/*.o misc/*.o surface/*.o volume/*.o
	- rm -f *.obj elements/*.obj gpath/*.obj grid/*.obj luawrap/*.obj misc/*.obj surface/*.obj volume/*.obj
	- rm -f $(DEMO_PROGRAMS)
	- rm -f $(TEST_PROGRAMS)
	- rm -f gpgrid*.vtk test*.vtk test*.stl test*.gz test*.svg test*.su2 test*.bin
	- rm -rf test_openFoam
	- cd $(LUA_DIR); make clean

all-clean:
	- make clean
	- rm -rf foamMesh foamMesh_with_rev_string.d
	- rm -rf $(BUILD_DIR)

# Specific targets, for individual artefacts.

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

foamMesh: foamMesh.d $(LIBLUA) $(GEOM_FILES) $(GEOM_LUA_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' foamMesh.d > foamMesh_with_rev_string.d

	$(DMD) $(DFLAGS_LUA) foamMesh_with_rev_string.d $(OF)foamMesh \
	$(GEOM_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

# The test artefacts.

vector3_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

vector3_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=vector3_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

projection_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

projection_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=projection_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

properties_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

properties_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=properties_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

line_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

line_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=line_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

arc_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

arc_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=arc_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

bezier_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

bezier_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=bezier_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

helix_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

helix_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=helix_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

nurbs_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

polyline_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

polyline_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=polyline_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

xspline_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

xsplinelsq_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

svgpath_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

coonspatch_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

coonspatch_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=coonspatch_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

aopatch_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

aopatch_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=aopatch_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

gmopatch_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

gmopatch_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=gmopatch_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

bezierpatch_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

bezierpatch_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=bezierpatch_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

controlpointpatch_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

nurbssurface_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

tfivolume_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

tfivolume_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=tfivolume_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

sweptsurfacevolume_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

sweptsurfacevolume_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=sweptsurfacevolume_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

twosurfacevolume_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

twosurfacevolume_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=twosurfacevolume_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

wedgevolume_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

wedgevolume_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=wedgevolume_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

slabvolume_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

slabvolume_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=slabvolume_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

nurbsvolume_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

sgrid_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

sgrid_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=sgrid_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

usgrid_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

usgrid_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=usgrid_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

univariatefunctions_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

svg_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

beziertrianglepatch_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

spherepatch_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

spherepatch_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=spherepatch_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

cubepatch_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

cubepatch_complex_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=cubepatch_test -version=complex_numbers \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

nurbs_utils_test: $(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -dip1008 -debug -g -w -version=$@ \
		$(GEOM_D_FILES) $(GZIP_FILES) $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)


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
