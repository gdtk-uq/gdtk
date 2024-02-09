GEOM_DIR ?= .

GEOM_D_FILES := \
	$(GEOM_DIR)/elements/nomenclature.d \
	$(GEOM_DIR)/elements/package.d \
	$(GEOM_DIR)/elements/projection.d \
	$(GEOM_DIR)/elements/properties.d \
	$(GEOM_DIR)/elements/vector3.d \
	\
	$(GEOM_DIR)/gpath/arc.d \
	$(GEOM_DIR)/gpath/bezier.d \
	$(GEOM_DIR)/gpath/gpath_utils.d \
	$(GEOM_DIR)/gpath/helix.d \
	$(GEOM_DIR)/gpath/line.d \
	$(GEOM_DIR)/gpath/modifiedpath.d \
	$(GEOM_DIR)/gpath/nurbs.d \
	$(GEOM_DIR)/gpath/package.d \
	$(GEOM_DIR)/gpath/path.d \
	$(GEOM_DIR)/gpath/polyline.d \
	$(GEOM_DIR)/gpath/polynomial.d \
	$(GEOM_DIR)/gpath/svgpath.d \
	$(GEOM_DIR)/gpath/xspline.d \
	$(GEOM_DIR)/gpath/xsplinelsq.d \
	\
	$(GEOM_DIR)/grid/grid.d \
	$(GEOM_DIR)/grid/package.d \
	$(GEOM_DIR)/grid/paver.d \
	$(GEOM_DIR)/grid/paver2d.d \
	$(GEOM_DIR)/grid/sgrid.d \
	$(GEOM_DIR)/grid/usgrid.d \
	\
	$(GEOM_DIR)/misc/kdtree.d \
	$(GEOM_DIR)/misc/nurbs_utils.d \
	$(GEOM_DIR)/misc/package.d \
	$(GEOM_DIR)/misc/sketch.d \
	$(GEOM_DIR)/misc/svg.d \
	$(GEOM_DIR)/misc/univariatefunctions.d \
	\
	$(GEOM_DIR)/surface/aopatch.d \
	$(GEOM_DIR)/surface/bezierpatch.d \
	$(GEOM_DIR)/surface/beziertrianglepatch.d \
	$(GEOM_DIR)/surface/channelpatch.d \
	$(GEOM_DIR)/surface/controlpointpatch.d \
	$(GEOM_DIR)/surface/coonspatch.d \
	$(GEOM_DIR)/surface/cubepatch.d \
	$(GEOM_DIR)/surface/gmopatch.d \
	$(GEOM_DIR)/surface/meshpatch.d \
	$(GEOM_DIR)/surface/nozzleexpansionpatch.d \
	$(GEOM_DIR)/surface/nurbssurface.d \
	$(GEOM_DIR)/surface/package.d \
	$(GEOM_DIR)/surface/parametricsurface.d \
	$(GEOM_DIR)/surface/ruledsurface.d \
	$(GEOM_DIR)/surface/spherepatch.d \
	$(GEOM_DIR)/surface/subrangedsurface.d \
	$(GEOM_DIR)/surface/sweptpathpatch.d \
	\
	$(GEOM_DIR)/volume/meshvolume.d \
	$(GEOM_DIR)/volume/nurbsvolume.d \
	$(GEOM_DIR)/volume/package.d \
	$(GEOM_DIR)/volume/parametricvolume.d \
	$(GEOM_DIR)/volume/slabvolume.d \
	$(GEOM_DIR)/volume/subrangedvolume.d \
	$(GEOM_DIR)/volume/sweptsurfacevolume.d \
	$(GEOM_DIR)/volume/tfivolume.d \
	$(GEOM_DIR)/volume/twosurfacevolume.d \
	$(GEOM_DIR)/volume/wedgevolume.d \
	\
	$(GEOM_DIR)/geometry_exception.d \
	\
	$(GEOM_DIR)/package.d \

# Look for the static library libplot
LIBPLOT := $(strip $(wildcard /usr/lib/libplot.a) \
                   $(wildcard $(LIBRARY_PATH)/libplot.a))
LIBPLOT_VERSION_STR :=
ifeq ($(findstring libplot,$(LIBPLOT)), libplot)
    $(warning Found libplot:$(LIBPLOT).)
    LIBPLOT_VERSION_STR := with_libplot
    GEOM_D_FILES := $(GEOM_D_FILES) $(GEOM_DIR)/misc/libplot.d
endif


GEOM_LUAWRAP_FILES := \
	$(GEOM_DIR)/luawrap/package.d \
	$(GEOM_DIR)/luawrap/luaunifunction.d \
	$(GEOM_DIR)/luawrap/luageom.d \
	$(GEOM_DIR)/luawrap/luanomenclature.d \
	$(GEOM_DIR)/luawrap/luagpath.d \
	$(GEOM_DIR)/luawrap/luagpath_utils.d \
	$(GEOM_DIR)/luawrap/luasurface.d \
	$(GEOM_DIR)/luawrap/luavolume.d \
	$(GEOM_DIR)/luawrap/luagrid.d \
	$(GEOM_DIR)/luawrap/luasgrid.d \
	$(GEOM_DIR)/luawrap/luausgrid.d \
	$(GEOM_DIR)/luawrap/luasketch.d

GEOM_FILES := $(GEOM_D_FILES) $(GEOM_LUAWRAP_FILES)

GEOM_LUA_FILES := $(GEOM_DIR)/foam-mesh.lua

