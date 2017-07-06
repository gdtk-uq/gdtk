# makefile for the grid utilities area

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd
PLATFORM ?= linux
INSTALL_DIR ?= $(HOME)/dgdinst
BUILD_DIR := ../../build

PROGRAMS := ugrid_partition

clean:
	- rm *.o
	- rm $(PROGRAMS)
	- rm $(BUILD_DIR)/bin/$(PROGRAMS)

install: $(PROGRAMS)
	cp -r $(BUILD_DIR)/bin/$(PROGRAMS) $(INSTALL_DIR)/bin

ugrid_partition: ugrid_partition.d
	$(DMD) -w -g ugrid_partition.d			
	cp $(PROGRAMS) $(BUILD_DIR)/bin/	

