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

ifeq ($(DMD), dmd)
    DFLAGS := -w -g -O -release -inline -boundscheck=off
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -g -O -release -enable-inlining -boundscheck=off -ffast-math -dip1008
endif
ifeq ($(DMD), ldc2)
    DFLAGS := -w -g -O -release -enable-inlining -boundscheck=off -ffast-math -dip1008
endif
ifeq ($(PLATFORM), macosx)
    DLINKFLAGS += -L-ld_classic
endif
PROGRAMS := ugrid_partition

ugrid_partition: ugrid_partition.d
	$(DMD) $(DFLAGS) ugrid_partition.d $(DLINKFLAGS)

clean:
	- rm -f *.o
	- rm -f $(PROGRAMS)
	- rm -f $(BUILD_DIR)/bin/$(PROGRAMS)

install: $(PROGRAMS)
	- mkdir -p $(BUILD_DIR)/bin
	cp $(PROGRAMS) $(BUILD_DIR)/bin/
	cp -r $(BUILD_DIR)/bin/$(PROGRAMS) $(INSTALL_DIR)/bin

