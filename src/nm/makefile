# makefile for the Numerical Methods package
#
# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= ldmd2

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

NTYPES_DIR := ../ntypes
NTYPES_FILES := $(NTYPES_DIR)/complex.d

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

include nm_files.mk

DEMO_PROGRAMS := ridder_demo \
	linesearch_demo \
	nelmin_demo nelmin_complex_demo \
	rungekutta_demo \
	newtoncotes_demo \
	gaussquad_demo \
	bbla_demo

LUA := ../../extern/lua-5.4.3
LIBLUA := $(LUA)/install/lib/liblua.a
LIBLUAPATH := $(LUA)/lib

ifeq ($(DMD), dmd)
    # DFLAGS := -w
    DFLAGS := -w -O -release -inline -boundscheck=off
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -O -release -inline -boundscheck=off
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif

# DIP1008 allows throwing of exceptions in @nogc code.
# See notes in src/eilmer/makefile.
DFLAGS += -dip1008 -preview=in

# ----------------------------------------------------------------------
# Here begins the list of targets, starting with the top-level actions.
#
# The default target is test, because that is most likely your interest
# if you are invoking this makefile from within the gas models directory.
# ----------------------------------------------------------------------

test: test_runner
	- exec ./test_runner --module nm

test-complex: test_runner_complex
	- exec ./test_runner_complex --module nm

test_runner: $(NM_FILES) $(NM_LUA_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA) ../util/test_runner.d
	@echo "Compiling test runner..."
	$(DMD) -of$@ -debug -g -w -dip1008 --unittest $(DLINKFLAGS) \
		../util/test_runner.d $(filter-out %/test_runner.d,$^)

test_runner_complex: $(NM_FILES) $(NM_LUA_FILES) $(NTYPES_FILES) $(UTIL_FILES) $(LIBLUA) ../util/test_runner.d
	@echo "Compiling test runner with complex numbers..."
	$(DMD) -of$@ -debug -g -w -dip1008 --unittest -version=complex_numbers $(DLINKFLAGS) \
		../util/test_runner.d $(filter-out %/test_runner.d,$^)

demo: ${DEMO_PROGRAMS}
	@echo "Demo codes built."

clean:
	- rm -f *.o *.obj
	- rm -f ${DEMO_PROGRAMS}
	- rm -f test_runner test_runner_complex
	- cd $(LUA); make clean

# Specific targets, for individual artefacts.

$(LIBLUA):
	cd $(LUA); make $(PLATFORM) local PLATFORM=$(PLATFORM)

# Some demo targets. These are left over from the early days of development.
# They should continue to work and provide a minimal example of how to get
# going with particular modules.

ridder_demo: ridder_demo.d ridder.d bracketing.d $(NTYPES_FILES)
	$(DMD) -g ridder_demo.d ridder.d bracketing.d $(NTYPES_FILES) $(DLINKFLAGS)

linesearch_demo: linesearch_demo.d linesearch.d $(NTYPES_FILES)
	$(DMD) -g linesearch_demo.d linesearch.d $(NTYPES_FILES) $(DLINKFLAGS)

nelmin_demo: nelmin_demo.d nelmin.d $(NTYPES_FILES) number.d
	$(DMD) -g nelmin_demo.d nelmin.d $(NTYPES_FILES) number.d $(DLINKFLAGS)

nelmin_complex_demo: nelmin_demo.d nelmin.d $(NTYPES_FILES) number.d
	$(DMD) -of$@ -g nelmin_demo.d nelmin.d $(NTYPES_FILES) number.d -version=complex_numbers $(DLINKFLAGS)

rungekutta_demo: rungekutta_demo.d rungekutta.d $(NTYPES_FILES)
	$(DMD) -g rungekutta_demo.d rungekutta.d $(NTYPES_FILES) $(DLINKFLAGS)

newtoncotes_demo: newtoncotes_demo.d newtoncotes.d $(NTYPES_FILES)
	$(DMD) -g newtoncotes_demo.d newtoncotes.d $(NTYPES_FILES) $(DLINKFLAGS)

gaussquad_demo: gaussquad_demo.d gaussquad.d $(NTYPES_FILES)
	$(DMD) -g gaussquad_demo.d gaussquad.d $(NTYPES_FILES) $(DLINKFLAGS)

bbla_demo: bbla_demo.d bbla.d $(NTYPES_FILES)
	$(DMD) -g -dip1008 bbla_demo.d bbla.d $(NTYPES_FILES) $(DLINKFLAGS)
