# makefile for the Numerical Methods package
#
# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= ldmd2
FLAVOUR ?= debug

DEMO_PROGRAMS := ridder_demo \
	linesearch_demo \
	nelmin_demo \
	rungekutta_demo \
	newtoncotes_demo \
	gaussquad_demo \
	bbla_demo

LUA := ../../extern/lua-5.4.3
LIBLUA := $(LUA)/install/lib/liblua.a

include ../universal.mk
include $(wildcard .reqs-*.mk)

d_previews := in dip1008
d_versions :=

imports := nm ntypes util
import_paths := ..

extra_deps := $(LIBLUA)

# ----------------------------------------------------------------------
# Here begins the list of targets, starting with the top-level actions.
#
# The default target is test, because that is most likely your interest
# if you are invoking this makefile from within the gas models directory.
# ----------------------------------------------------------------------

test-all: test-real test-complex

test-real: test_runner_real
	- exec ./test_runner_real --module nm

test-complex: test_runner_complex
	- exec ./test_runner_complex --module nm

test_runner_%: FLAVOUR = debug
test_runner_%: compiler_flags += -unittest

test_runner_real: ../util/test_runner.d package.d $(extra_deps)
	@echo "Compiling test runner..."
	$(DMD) -of=$@ -makedeps=.reqs-$@.mk $(compiler_flags) $< package.d $(extra_deps)

test_runner_complex: d_versions += complex_numbers
test_runner_complex: ../util/test_runner.d package.d $(extra_deps)
	@echo "Compiling test runner with complex numbers..."
	$(DMD) -of=$@ -makedeps=.reqs-$@.mk $(compiler_flags) $< package.d $(extra_deps)

demo: $(DEMO_PROGRAMS) nelmin_complex_demo
	@echo "Demo codes built."

# Template for demo DEMO_PROGRAMS
$(DEMO_PROGRAMS): %: %.d
	$(DMD) -of=$@ -makedeps=.reqs-$@.mk $(compiler_flags) $<

# Specialise the nelmin complex version slightly 
nelmin_complex_demo: d_versions += complex_numbers
nelmin_complex_demo: nelmin_demo.d
	$(DMD) -of=$@ -makedeps=.reqs-$@.mk $(compiler_flags) $<

clean:
	- rm -f *.o *.obj
	- rm -f .reqs-*.mk
	- rm -f $(DEMO_PROGRAMS) nelmin_complex_demo
	- rm -f test_runner test_runner_complex
	- cd $(LUA); make clean

# Specific targets, for individual artefacts.

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

$(LIBLUA):
	cd $(LUA); make $(PLATFORM) local PLATFORM=$(PLATFORM)

.PHONY: test-all test-real test-complex clean
