# makefile for the geometry package
# Builds all of the demonstration programs.

include geom_files.mk

DMD ?= ldc2

INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR := ../../build
REPO_DIR := $(shell cd ../../; pwd)
ifneq ("$(wildcard $(REPO_DIR)/.git)","")
    # Repository exists
    REVISION_STRING := $(shell git rev-parse --short HEAD)
else
    REVISION_STRING := unknown
endif

demo_files := $(wildcard demos/*)
demo_targets := $(basename $(notdir $(demo_files)))

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a

# ----------------------------------------------------------------------------
# The default target is test-all, because that is most likely your interest
# if you are invoking this makefile from within the geom models directory.
# ----------------------------------------------------------------------------

default: test-all

play:
	@echo "Found $(LIBPLOT)"

install: foamMesh grid2pato $(OPENFOAM_DIR)/foam-mesh.lua $(OPENFOAM_DIR)/openfoam.lua
	- mkdir -p $(INSTALL_DIR)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share
	cp foamMesh $(BUILD_DIR)/bin/
	cp grid2pato $(BUILD_DIR)/bin/
	cp $(OPENFOAM_DIR)/foam-mesh.lua $(BUILD_DIR)/bin/
	cp $(OPENFOAM_DIR)/grid2pato.lua $(BUILD_DIR)/bin/
	cp $(OPENFOAM_DIR)/openfoam.lua $(BUILD_DIR)/lib/
	cp -r $(OPENFOAM_DIR)/foamMesh-templates $(BUILD_DIR)/share/
	cp -r $(OPENFOAM_DIR)/pato-templates $(BUILD_DIR)/share/
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

clean:
	- rm -f *.o **/*.o
	- rm -f *.obj **/*.obj
	- rm -f *.deps **/*.deps
	- rm -rf foamMesh grid2pato
	- rm -f $(demo_targets)
	- rm -f test_runner_real test_runner_complex
	- rm -f gpgrid*.vtk test*.vtk test*.stl test*.gz test*.svg test*.su2 test*.bin
	- rm -rf test_openFoam
	- rm **/*_with_rev_string.d
	- rm $(BUILD_DIR)/lib/openfoam.lua
	- $(MAKE) -C $(LUA_DIR) clean

all-clean: clean
	- rm -rf foamMesh
	- rm -rf grid2pato
	- rm -rf $(BUILD_DIR)

include ../common.mk

d_previews := in dip1008
d_versions := flavour_$(FLAVOUR)

# External dependencies
imports := gzip
import_paths := ../extern/gzip

imports += geom nm ntypes util
import_paths += ..

extra_deps := $(LIBLUA)

ifeq ($(findstring libplot,$(LIBPLOT_VERSION_STR)), libplot)
	d_versions += $(LIBPLOT_VERSION_STR)
	linker_libs += plot
endif

# ----------------------------------------------------------------------
# Here begins the list of targets, starting with the top-level actions.
# ----------------------------------------------------------------------

include $(wildcard *.deps)

.PHONY: test-all test-real test-complex
test-all: test-real test-complex

test-real: test_runner_real
	@echo "Running real tests"
	- exec ./test_runner_real --module geom

test-complex: test_runner_complex
	@echo "Running complex tests"
	- exec ./test_runner_complex --module geom

test_runner_%: FLAVOUR = debug
test_runner_%: compiler_flags += -unittest
test_runner_complex: d_versions += complex_numbers
test_runner_complex test_runner_real: ../util/test_runner.d package.d
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

# Specific targets, for individual artefacts.

%_with_rev_string.d: %.d
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' $< > $@

foamMesh grid2pato: %: $(OPENFOAM_DIR)/%_with_rev_string.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)


# The original demo artefacts.

demo: $(demo_targets)
	@echo "Demo codes built."

$(demo_targets): %: demos/%.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

$(LIBLUA):
	$(MAKE) -C $(LUA_DIR) guess local
