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
LIBLUAPATH := $(LUA_DIR)/install/lib
LIBLUA := $(LIBLUAPATH)/liblua.a

EQC_DIR := ../extern/eqc/source
LIBEQC := $(EQC_DIR)/libeqc.a

include ../universal.mk
include $(wildcard *.deps)

d_previews := in dip1008
d_versions :=

# External dependencies
imports := eqc
import_paths := ../extern/eqc/source

imports += gasdyn nm ntypes util
import_paths += ..

extra_deps := $(LIBLUA)

.SECONDEXPANSION:

# Specific requirements for gasflow demos
gasflow_demo luagasflow_demo: imports += gas kinetics
gasflow_demo luagasflow_demo: extra_deps += $(LIBEQC)

# Template for all demo targets
# As the demo targetNames are the same as their source-file,
# these targets are *super* simple
$(demo_targets): %: %.d $$(extra_deps)
	$(DMD) -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

GAS_DIR := ../gas
KINETICS_DIR := ../kinetics

clean:
	- rm *.o *.obj
	- rm *.deps
	- rm $(demo_targets)
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(EQC_DIR); make clean

$(LIBLUA):
	cd $(LUA_DIR); make guess local

$(LIBEQC):
	cd $(EQC_DIR); make
