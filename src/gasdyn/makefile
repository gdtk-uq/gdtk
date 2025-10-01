# makefile for gasdyn
# Builds all of the demonstration programs.
#
# Peter J. September 2015
#          October 2016 bring into line with makefile from geom
#

demo_targets := idealgasflow_demo luaidealgasflow_demo \
	gasflow_demo luagasflow_demo

demo: $(demo_targets)
	@echo "Demo codes built."

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a

include $(wildcard .reqs-*.mk)
include ../universal.mk

d_previews := in dip1008
d_versions :=

# External dependencies
imports := eqc
import_paths := ../extern/eqc/source

imports += gasdyn nm ntypes util
import_paths += ..

extra_deps := $(LIBLUA)

# Specific requirements for gasflow demos
gasflow_demo luagasflow_demo: $(LIBEQC)
gasflow_demo luagasflow_demo: imports += gas kinetics
gasflow_demo luagasflow_demo: extra_deps += $(LIBEQC)

# Template for all demo targets
$(demo_targets): %: %.d $(extra_deps)
	$(DMD) -of=$@ -makedeps=.reqs-$@.mk $(compiler_flags) $< $(extra_deps)

GAS_DIR := ../gas
KINETICS_DIR := ../kinetics

clean:
	- rm *.o *.obj
	- rm .reqs-*.mk
	- rm $(demo_targets)
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(EQC_DIR); make clean

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

$(LIBEQC):
	cd $(EQC_DIR); make
