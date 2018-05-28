# makefile for the Numerical Methods package
#
# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd
PLATFORM ?= linux

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

include nm_files.mk

DEMO_PROGRAMS := ridder_demo \
	linesearch_demo \
	nelmin_demo \
	rungekutta_demo \
	newtoncotes_demo \
	gaussquad_demo \
	bbla_demo 

TEST_PROGRAMS := bbla_test \
	bracketing_test bracketing_complex_test \
	gaussquad_test gaussquad_complex_test \
	linesearch_test \
	newtoncotes_test \
	ridder_test \
	brent_test brent_complex_test \
	rungekutta_test rungekutta_complex_test \
	smla_test \
	rsla_test \
	luabbla_test

LUA := ../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib
ifeq ($(DMD), dmd)
    # DFLAGS := -w
    DFLAGS := -w -O -release -inline -boundscheck=off
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -O -release -inline -boundscheck=off
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), gdc)
    DFLAGS := -Wall -O2 -frelease -fno-debug
    OF := -o
    DLINKFLAGS := -L$(LIBLUAPATH) $(LIBLUA) -ldl
endif

# ----------------------------------------------------------------------
# Here begins the list of targets, starting with the top-level actions.
#
# The default target is test, because that is most likely your interest
# if you are invoking this makefile from within the gas models directory.
# ----------------------------------------------------------------------

test: ${TEST_PROGRAMS}
	tclsh nm-package-test.tcl

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}
	- rm ${TEST_PROGRAMS}
	- cd $(LUA); make clean

# Specific targets, for individual artefacts.

$(LIBLUA):
	cd $(LUA); make linux local

bbla_test: bbla.d $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=$@ \
		bbla.d $(UTIL_DIR)/msg_service.d

bracketing_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

bracketing_complex_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=bracketing_test -version=complex_numbers \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

gaussquad_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

gaussquad_complex_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=gaussquad_test -version=complex_numbers \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

linesearch_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

newtoncotes_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

ridder_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

rungekutta_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

rungekutta_complex_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=rungekutta_test -version=complex_numbers \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

brent_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

brent_complex_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=brent_test -version=complex_numbers \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

smla_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

rsla_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -g -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

luabbla_test: luabbla.d bbla.d $(UTIL_FILES) $(LIBLUA) $(UTIL_FILES)
	dmd -of$@ -debug -g -w -version=$@ \
		luabbla.d bbla.d $(UTIL_FILES) $(DLINKFLAGS)

# Some demo targets. These are left over from the early days of development.
# They should continue to work and provide a minimal example of how to get
# going with particular modules.

ridder_demo: ridder_demo.d ridder.d bracketing.d
	dmd -g ridder_demo.d ridder.d bracketing.d

linesearch_demo: linesearch_demo.d linesearch.d
	dmd -g linesearch_demo.d linesearch.d

nelmin_demo: nelmin_demo.d nelmin.d
	dmd -g nelmin_demo.d nelmin.d

rungekutta_demo: rungekutta_demo.d rungekutta.d
	dmd -g rungekutta_demo.d rungekutta.d

newtoncotes_demo: newtoncotes_demo.d newtoncotes.d
	dmd -g newtoncotes_demo.d newtoncotes.d

gaussquad_demo: gaussquad_demo.d gaussquad.d
	dmd -g gaussquad_demo.d gaussquad.d

bbla_demo: bbla_demo.d bbla.d
	dmd -g bbla_demo.d bbla.d
