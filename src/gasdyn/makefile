# makefile for gasdyn
#
# Peter J. September 2015
#
# Adapted from the Eilmer4 makefile so it has a bit of stuff left in for future use.  
# (For example, there's the Lua build, for when we connect up the fancier gas models.)

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd

PROGRAMS := 

DEMO_PROGRAMS := idealgasflow_demo

NM_DIR := ../nm
NM_FILES := $(NM_DIR)/linesearch.d

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
    #     DFLAGS := -debug -gc -w -unittest
    # (2) For profiling the debugged code.
    #     DFLAGS := -profile -vgc -w -O -release -inline -boundscheck=off
    # (3) For releasing a fast code.
    #     DFLAGS := -w -O -release -inline -boundscheck=off
    DFLAGS := -w -O -debug -gc -unittest \
	-I.. -I$(NM_DIR) -I$(GEOM_DIR) -I$(GZIP_DIR)
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -O -release -inline -boundscheck=off \
	-I.. -I$(NM_DIR) -I$(GEOM_DIR) -I$(GZIP_DIR)
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), gdc)
    DFLAGS := -Wall -O2 -frelease -fno-debug \
	-I.. -I$(NM_DIR) -I$(GEOM_DIR) -I$(GZIP_DIR)
    OF := -o
    DLINKFLAGS := -L$(LIBLUAPATH) $(LIBLUA) -ldl
endif

default: $(PROGRAMS) $(DEMO_PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "gasdyn demo programs built."

install: $(PROGRAMS)
	#@echo "Installing to $(INSTALL_DIR)"
	#cp $(PROGRAMS) $(INSTALL_DIR)
	@echo "Do nothing"

demo: $(DEMO_PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Demo codes built."

clean:
	- rm *.o
	- rm $(DEMO_PROGRAMS) $(PROGRAMS)
	- rm trace.log trace.def
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a

$(LIBLUA):
	cd $(LUA_DIR); make linux local

idealgasflow_demo: idealgasflow_demo.d idealgasflow.d $(NM_FILES) $(LIBLUA)
	$(DMD) $(DFLAGS) $(OF)idealgasflow_demo \
		idealgasflow_demo.d idealgasflow.d $(NM_FILES) $(DLINKFLAGS)
