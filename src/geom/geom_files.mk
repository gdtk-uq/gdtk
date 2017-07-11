GEOM_DIR ?= .
GEOM_FILES := $(GEOM_DIR)/geom.d \
	$(GEOM_DIR)/univariatefunctions.d \
	$(GEOM_DIR)/luaunifunction.d \
	$(GEOM_DIR)/gpath.d \
	$(GEOM_DIR)/surface.d \
	$(GEOM_DIR)/volume.d \
	$(GEOM_DIR)/luageom.d \
	$(GEOM_DIR)/luagpath.d \
	$(GEOM_DIR)/luasurface.d \
	$(GEOM_DIR)/luavolume.d \
	$(GEOM_DIR)/grid.d \
	$(GEOM_DIR)/luagrid.d \
	$(GEOM_DIR)/sgrid.d \
	$(GEOM_DIR)/luasgrid.d \
	$(GEOM_DIR)/usgrid.d \
	$(GEOM_DIR)/luausgrid.d \
	$(GEOM_DIR)/paver.d \
	$(GEOM_DIR)/svg.d \
	$(GEOM_DIR)/sketch.d \
	$(GEOM_DIR)/luasketch.d

GEOM_LUA_FILES := $(GEOM_DIR)/foam-mesh.lua

