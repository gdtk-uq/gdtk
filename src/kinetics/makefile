# makefile for kinetics module
# Can be used to perform unit tests
# and build stand-alone programs.
PLATFORM ?= linux
WITH_DVODE ?= 0
INSTALL_DIR ?= $(HOME)/e3bin
BUILD_DIR ?= ../../build

TEST_PROGRAMS := chemistry_update_test \
		equilibrium_update_test \
		rate_constant_test \
		reaction_test \
		reaction_mechanism_test \
		two_temperature_argon_kinetics_test \
		two_temperature_argon_with_ideal_gas_test \
		electronically_specific_kinetics_test \
		two_temperature_gasgiant_kinetics_test

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk

CEQ_DIR := ../extern/ceq/source
LIBCEQ := $(CEQ_DIR)/libceq.a
include $(CEQ_DIR)/ceq_files.mk

include kinetics_files.mk

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := ${LUA_DIR}/lib/liblua.a
LIBLUAPATH := ${LUA_DIR}/lib

DLINKFLAGS :=  -w -L-L${LIBLUAPATH} -L-llua -L-ldl

ifeq ($(WITH_DVODE),1)
    DFLAGS += -version=with_dvode
    DLINKFLAGS += -L-lgfortran
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

test: $(TEST_PROGRAMS)
	tclsh kinetics-package-test.tcl

clean:
	- rm -f *.o *.mod
	- rm -f $(TEST_PROGRAMS)
	- rm -f two_temperature_argon_kinetics_test_results.data
	- rm -rf $(BUILD_DIR)/*
	- cd $(LUA_DIR); make clean
	- cd $(CEQ_DIR); make clean
	- rm -f $(GAS_DIR)/ideal_air_fortran.mod $(GAS_DIR)/ideal_air_fortran.o

${LIBLUA}:
	cd ${LUA_DIR}; make $(PLATFORM) local

$(LIBCEQ):
	cd $(CEQ_DIR); make

chemistry_update_test: $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	dmd -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

equilibrium_update_test: $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	dmd -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

rate_constant_test: $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	dmd -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

reaction_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	dmd -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

reaction_mechanism_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	dmd -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

two_temperature_argon_kinetics_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	dmd -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

two_temperature_argon_with_ideal_gas_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	dmd -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

pseudo_species_kinetics_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	dmd -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

electronically_specific_kinetics_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	dmd -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)

two_temperature_gasgiant_kinetics_test:  $(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(LIBCEQ) $(LIBLUA)
	dmd -of$@ -debug -g -dip1008 -version=$@ \
		$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) \
		$(LIBCEQ) $(LIBLUA) $(DLINKFLAGS)


#--------------------------------------------------------------------------
# Fortran bits

$(GAS_DIR)/ideal_air_fortran.o: $(GAS_DIR)/ideal_air_fortran.f
	cd $(GAS_DIR); gfortran -c -ffree-form ideal_air_fortran.f

pseudosp_solve_ode.o: pseudosp_solve_ode.f90
	gfortran -c pseudosp_solve_ode.f90

pseudosp_rhs_ode.o: $(KINETICS_DIR)/pseudosp_rhs_ode.f90
	gfortran -c pseudosp_rhs_ode.f90

dvode.o: $(KINETICS_DIR)/dvode.f
	gfortran -c dvode.f
