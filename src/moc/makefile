# makefile for MOC
#
# Peter J. and Momar Hughes, September 2015
#
# Adapted from the Eilmer4 makefile so it has a bit of stuff left in
# for future use.  (For example, there's the Lua build.)

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd

PROGRAMS := moc

DEMO_PROGRAMS := kernel_demo \
		wall_demo \
		anderson_example_11_1

MOC_DIR := .
MOC_FILES := $(MOC_DIR)/kernel.d \
	$(MOC_DIR)/unitproc.d \
	$(MOC_DIR)/wall.d

NM_DIR := ../nm
NM_FILES := $(NM_DIR)/bbla.d \
	$(NM_DIR)/linesearch.d \
	$(NM_DIR)/ridder.d

GEOM_DIR := ../geom
GEOM_FILES := $(GEOM_DIR)/geom.d \
	$(GEOM_DIR)/gpath.d

GASDYN_DIR := ../gasdyn
GASDYN_FILES := $(GASDYN_DIR)/idealgasflow.d

GAS_DIR := ../gas
GAS_FILES := $(GAS_DIR)/gas_model.d \
	$(GAS_DIR)/physical_constants.d

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/msg_service.d

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := $(LUA_DIR)/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib

# The install destination.
INSTALL_DIR ?= $(HOME)/dgdinst

REVISION_STRING := $(shell hg identify --id --num --branch --tags)

ifeq ($(DMD), dmd)
    # (1) For debug build.
    #     DFLAGS := -debug -g -w -unittest
    # (2) For profiling the debugged code.
    #     DFLAGS := -profile -vgc -w -O -release -inline -boundscheck=off
    # (3) For releasing a fast code.
    #     DFLAGS := -w -O -release -inline -boundscheck=off
    DFLAGS := -w -O -debug -g -unittest \
	-I.. -I$(NM_DIR) -I$(GEOM_DIR) -I$(GAS_DIR) -I$(UTIL_DIR) -I$(GZIP_DIR)
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -O -release -inline -boundscheck=off \
	-I.. -I$(NM_DIR) -I$(GEOM_DIR) -I$(GAS_DIR) -I$(UTIL_DIR) -I$(GZIP_DIR)
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), gdc)
    DFLAGS := -Wall -O2 -frelease -fno-debug \
	-I.. -I$(NM_DIR) -I$(GEOM_DIR) -I$(GAS_DIR) -I$(UTIL_DIR) -I$(GZIP_DIR)
    OF := -o
    DLINKFLAGS := -L$(LIBLUAPATH) $(LIBLUA) -ldl
endif

default: $(PROGRAMS) $(DEMO_PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "moc demo programs built."
	@echo "One day, this will also build the main moc program."

install: $(PROGRAMS)
	@echo "Installing to $(INSTALL_DIR)"
	cp $(PROGRAMS) $(INSTALL_DIR)

demo: $(DEMO_PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Demo codes built."

clean:
	- rm *.o
	- rm $(DEMO_PROGRAMS) $(PROGRAMS)
	- rm trace.log trace.def
	- rm anderson_11_1.txt kernel_demo.txt moc_wall_demo.txt
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a

$(LIBLUA):
	cd $(LUA_DIR); make linux local

moc: 
	@echo "Place-holder for the main moc program."

kernel_demo: kernel_demo.d kernel.d \
	$(GEOM_FILES) $(GASDYN_FILES) $(GAS_FILES) \
	$(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(DFLAGS) $(OF)kernel_demo \
		kernel_demo.d kernel.d \
		$(GEOM_FILES) $(GASDYN_FILES) $(GAS_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

wall_demo: wall_demo.d wall.d kernel.d\
	$(GEOM_FILES) $(GASDYN_FILES) $(GAS_FILES) \
	$(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(DFLAGS) $(OF)wall_demo \
		wall_demo.d wall.d kernel.d\
		$(GEOM_FILES) $(GASDYN_FILES) $(GAS_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

anderson_example_11_1: anderson_example_11_1.d $(MOC_FILES)\
	$(GEOM_FILES) $(GASDYN_FILES) $(GAS_FILES) \
	$(UTIL_FILES) $(NM_FILES) $(LIBLUA)
	$(DMD) $(DFLAGS) $(OF)anderson_example_11_1 \
		anderson_example_11_1.d $(MOC_FILES) \
		$(GEOM_FILES) $(GASDYN_FILES) $(GAS_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)  

