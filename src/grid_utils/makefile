# makefile for the grid utilities area

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= ldc2

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR ?= ../../build

FLAVOUR ?= fast
PROGRAMS := ugrid_partition

include ../universal.mk # Generates compiler_flags for us
include $(wildcard .reqs-*.mk)

# Template for all d-programs
$(PROGRAMS): %: %.d
	$(DMD) -of=$@ -makedeps=.reqs-$@.mk $(compiler_flags) $<

clean:
	- rm -f *.o .reqs-*.mk
	- rm -f $(PROGRAMS)
	- rm -f $(BUILD_DIR)/bin/$(PROGRAMS)

install: $(PROGRAMS)
	- mkdir -p $(BUILD_DIR)/bin
	cp $(PROGRAMS) $(BUILD_DIR)/bin/
	cp -r $(BUILD_DIR)/bin/$(PROGRAMS) $(INSTALL_DIR)/bin

