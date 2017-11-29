GEOM_DIR ?= .

GEOM_D_FILES := \
	$(GEOM_DIR)/package.d \
	$(GEOM_DIR)/elements/package.d \
	$(GEOM_DIR)/elements/nomenclature.d \
	$(GEOM_DIR)/elements/vector3.d \
	$(GEOM_DIR)/elements/properties.d \
	$(GEOM_DIR)/elements/projection.d \
	$(GEOM_DIR)/gpath/package.d \
	$(GEOM_DIR)/gpath/path.d \
	$(GEOM_DIR)/gpath/line.d \
	$(GEOM_DIR)/gpath/arc.d \
	$(GEOM_DIR)/gpath/helix.d \
	$(GEOM_DIR)/gpath/bezier.d \
	$(GEOM_DIR)/gpath/polynomial.d \
	$(GEOM_DIR)/gpath/polyline.d \
	$(GEOM_DIR)/gpath/modifiedpath.d \
	$(GEOM_DIR)/surface/package.d \
	$(GEOM_DIR)/surface/surface.d \
	$(GEOM_DIR)/volume/package.d \
	$(GEOM_DIR)/volume/volume.d \
	$(GEOM_DIR)/grid/package.d \
	$(GEOM_DIR)/grid/grid.d \
	$(GEOM_DIR)/grid/sgrid.d \
	$(GEOM_DIR)/grid/usgrid.d \
	$(GEOM_DIR)/grid/paver.d \
	$(GEOM_DIR)/misc/package.d \
	$(GEOM_DIR)/misc/univariatefunctions.d \
	$(GEOM_DIR)/misc/svg.d \
	$(GEOM_DIR)/misc/sketch.d

GEOM_LUAWRAP_FILES := \
	$(GEOM_DIR)/luawrap/package.d \
	$(GEOM_DIR)/luawrap/luaunifunction.d \
	$(GEOM_DIR)/luawrap/luageom.d \
	$(GEOM_DIR)/luawrap/luagpath.d \
	$(GEOM_DIR)/luawrap/luasurface.d \
	$(GEOM_DIR)/luawrap/luavolume.d \
	$(GEOM_DIR)/luawrap/luagrid.d \
	$(GEOM_DIR)/luawrap/luasgrid.d \
	$(GEOM_DIR)/luawrap/luausgrid.d \
	$(GEOM_DIR)/luawrap/luasketch.d

GEOM_FILES := $(GEOM_D_FILES) $(GEOM_LUAWRAP_FILES)

GEOM_LUA_FILES := $(GEOM_DIR)/foam-mesh.lua

