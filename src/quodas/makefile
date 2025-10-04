# makefile: quodas
# Kyle A. Damm (27-01-2020)

DMD := ldc2
FLAVOUR ?= fast
WITH_COMPLEX_NUMBERS ?= 0

include ../common.mk
include $(wildcard *.deps)

d_versions := \
	$(if $(filter 1,$(WITH_COMPLEX_NUMBERS)), complex_numbers)

imports := nm ntypes util
import_paths := ..

# Because of the project structure, D compilers struggle to
# correctly import the quodas project. We feed in all
# the source files manually.
source_files = $(wildcard source/*.d)

quodas: $(source_files)
	- $(DMD) -of=$@ -makedeps=$@.deps $(compiler_flags) $^
	$(call cleanup-deps)

clean:
	- rm *.o
	- rm *.deps
	- rm quodas
