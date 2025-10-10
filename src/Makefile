INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR := ../build

LUA_DIR := ../extern/lua-5.4.3
LPEG_DIR := ../extern/lpeg-1.0.2

export FLAVOUR ?= debug

DIRS := \
	chicken \
	gas \
	gasdyn \
	geom \
	grid_utils \
	kinetics \
	l1d \
	lmr \
	nenzf1d \
	nm \
	ntypes \
	puffin \
	slf \
	util

clean:
	@for dir in $(DIRS); do \
		echo "Cleaning directory $$dir/"; \
		make -s -C $$dir clean; \
	done
	# For some reason geom thinks it's special
	make -s -C geom all-clean
	rm -rf $(BUILD_DIR)

uninstall:
	rm -rf $(INSTALL_DIR)
	
prep-install:
	- mkdir -p $(INSTALL_DIR)/bin
	- mkdir -p $(INSTALL_DIR)/etc
	- mkdir -p $(INSTALL_DIR)/lib
	- mkdir -p $(INSTALL_DIR)/share

# This should really be done at the program-build level
# Only artifacts that are generated should be moved here.
# Anything else should just be done directly in "install"
prep-build-dir:
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share

install-lua: prep-install
	make -C $(LUA_DIR) guess local
	cp $(LUA_DIR)/install/bin/* $(INSTALL_DIR)/bin

install-lpeg: prep-install
	make -C $(LPEG_DIR) guess
	cp $(LPEG_DIR)/lpeg.so $(INSTALL_DIR)/lib

install-gdtklib: prep-install
	# Just copy straight to install-dir?
	cp -r ./lib/* $(INSTALL_DIR)/lib

install-gas: prep-install install-lua
	make -C ./gas prep-gas libgas
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

install-geom: OPENFOAM_DIR := openfoam
install-geom: prep-install prep-build-dir install-lua
	make -C ./geom foamMesh grid2pato
	cp ./geom/foamMesh ./geom/grid2pato $(BUILD_DIR)/bin/
	cp ./geom/$(OPENFOAM_DIR)/foam-mesh.lua $(BUILD_DIR)/bin/
	cp ./geom/$(OPENFOAM_DIR)/grid2pato.lua $(BUILD_DIR)/bin/
	cp ./geom/$(OPENFOAM_DIR)/openfoam.lua $(BUILD_DIR)/lib/
	cp -r ./geom/$(OPENFOAM_DIR)/foamMesh-templates $(BUILD_DIR)/share/
	cp -r ./geom/$(OPENFOAM_DIR)/pato-templates $(BUILD_DIR)/share/
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

install-gridutils: PROGRAMS := ugrid_partition
install-gridutils: prep-install prep-build-dir
	make -C ./grid_utils $(PROGRAMS)
	cp $(addprefix ./grid_utils/,$(PROGRAMS)) $(BUILD_DIR)/bin/
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

install-kinetics: PROGRAMS := prep-chem prep-kinetics chemkin2eilmer
install-kinetics: prep-install prep-build-dir install-lua install-lpeg
	make -C ./kinetics $(PROGRAMS)
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)
	
install-l1d: PROGRAMS := l1d4 l1d4-prep l1d4-xt-viewer
install-l1d: prep-install prep-build-dir install-lua install-gdtklib install-gas install-kinetics
	make -C ./l1d $(PROGRAMS)
	cp $(addprefix ./l1d/,$(PROGRAMS)) $(BUILD_DIR)/bin
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

install-lmr: PROGRAMS := lmr lmr-debug
install-lmr: SUB_PROGRAMS := lmr-run lmrZ-run
install-lmr: PY_PROGRAMS := lmr-verify
install-lmr: AUX_PROGRAMS := lmr-complete gdtk-module
install-lmr: prep-install prep-build-dir install-lua install-lpeg install-gdtklib install-gas install-kinetics install-gridutils
	make -C ./lmr $(PROGRAMS) $(SUB_PROGRAMS) $(PY_PROGRAMS) $(AUX_PROGRAMS)
	@echo "Installing to $(INSTALL_DIR)"
	cp $(addprefix ./lmr/,$(PROGRAMS)) $(INSTALL_DIR)/bin
	cp $(addprefix ./lmr/,$(SUB_PROGRAMS)) $(INSTALL_DIR)/bin
	cp $(addprefix ./lmr/,$(PY_PROGRAMS)) $(INSTALL_DIR)/bin
	# Where is AUX_PROGRAMS?
	cp ./lmr/lua-modules/*.lua $(INSTALL_DIR)/lib/
	cp ./lmr/lmr.cfg $(INSTALL_DIR)/etc/
	cp ./lmr/gdtk-module $(INSTALL_DIR)/share
	cp ./lmr/lmr-complete.sh $(INSTALL_DIR)/share
	cp ./lmr/share/* $(INSTALL_DIR)/share
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

install-nenzf1d: PROGRAMS := nenzf1d
install-nenzf1d: prep-install prep-build-dir install-gdtklib install-gas install-kinetics
	make -C ./nenzf1d $(PROGRAMS) libnenzf1d.so
	cp $(addprefix ./nenzf1d/,$(PROGRAMS)) $(BUILD_DIR)/bin
	cp ./nenzf1d/nenzf1d.py $(BUILD_DIR)/lib
	cp ./nenzf1d/libnenzf1d.so $(BUILD_DIR)/lib
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

install-slf: PROGRAMS := slf libslf.so
install-slf: prep-install prep-build-dir install-gdtklib install-gas install-kinetics
	make -C ./slf $(PROGRAMS)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib/gdtk
	- cp ./slf/slf $(BUILD_DIR)/bin
	- cp ./slf/libslf.so $(BUILD_DIR)/lib
	- cp ./slf/slf.py $(BUILD_DIR)/lib/gdtk
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)
