# makefile for kinetics module
# Can be used to perform unit tests
# and build stand-alone programs.

DMD := ldmd2

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR ?= ../../build

TEST_PROGRAMS := chemistry_update_test \
		equilibrium_update_test \
		rate_constant_test \
		reaction_test \
		reaction_mechanism_test \
		two_temperature_argon_kinetics_test \
		vib_specific_co_kinetics_test \
		vib_specific_co_mixture_kinetics_test \
		two_temperature_argon_with_ideal_gas_test \
		two_temperature_gasgiant_kinetics_test

NTYPES_DIR := ../ntypes
NTYPES_FILES := $(NTYPES_DIR)/complex.d

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk
LIBGASF := $(GAS_DIR)/libgasf.a

CEQ_DIR := ../extern/ceq/source
LIBCEQ := $(CEQ_DIR)/libceq.a
include $(CEQ_DIR)/ceq_files.mk

include kinetics_files.mk

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := ${LUA_DIR}/install/lib/liblua.a
LIBLUAPATH := ${LUA_DIR}/lib

DLINKFLAGS :=  -w -L-ldl
ifeq ($(PLATFORM), macosx)
    DLINKFLAGS += -L-ld_classic
endif

build-prep-chem: prep_chem.lua reaction.lua lex_elems.lua
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	cp reaction.lua lex_elems.lua $(BUILD_DIR)/lib/
	cp prep_chem.lua $(BUILD_DIR)/bin/prep-chem; chmod +x $(BUILD_DIR)/bin/prep-chem

build-chemkin2eilmer: chemkin2eilmer.lua lex_elems.lua reaction.lua
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	cp reaction.lua lex_elems.lua $(BUILD_DIR)/lib/
	cp chemkin2eilmer.lua $(BUILD_DIR)/bin/chemkin2eilmer; chmod +x $(BUILD_DIR)/bin/chemkin2eilmer

build-prep-kinetics: prep_kinetics.lua mechanism.lua lex_elems.lua
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	cp mechanism.lua lex_elems.lua $(BUILD_DIR)/lib/
	cp prep_kinetics.lua $(BUILD_DIR)/bin/prep-kinetics; chmod +x $(BUILD_DIR)/bin/prep-kinetics

test: $(TEST_PROGRAMS)
	tclsh kinetics-package-test.tcl

clean:
	- rm -f *.o *.mod *.obj
	- rm -f $(TEST_PROGRAMS)
	- rm -f two_temperature_argon_kinetics_test_results.data
	- rm -rf $(BUILD_DIR)/*
	- cd $(LUA_DIR); make clean
	- cd $(CEQ_DIR); make clean

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

$(LIBCEQ):
	cd $(CEQ_DIR); make

chemistry_update_test: $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

equilibrium_update_test: $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

rate_constant_test: $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

reaction_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

reaction_mechanism_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

two_temperature_argon_kinetics_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

vib_specific_co_kinetics_test: $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

vib_specific_co_mixture_kinetics_test: $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

two_temperature_argon_with_ideal_gas_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

two_temperature_gasgiant_kinetics_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)
