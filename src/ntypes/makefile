# makefile for the Numerical Types package
#
# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd

include ../common.mk
include $(wildcard *.deps)

d_previews := in dip1008
d_versions := complex_numbers

imports := nm ntypes util
import_paths := ..

.PHONY: test
test: test_runner
	./test_runner --module ntypes

test_runner: FLAVOUR = debug
test_runner: compiler_flags += -unittest
test_runner: ../util/test_runner.d complex.d
	@echo "Compiling test runner..."
	$(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

clean:
	- rm -f *.o *.deps
	- rm -f test_runner
