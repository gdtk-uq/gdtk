# makefile for Lagrangian 1D flow solver.
# PJ 2020-04-05

DMD ?= ldc2

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug
WITH_FPE ?= 1

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

PROGRAMS := l1d4 l1d4-prep l1d4-xt-viewer

L1D_FILES := config.d simcore.d tube.d gasslug.d lcell.d \
	piston.d valve.d endcondition.d misc.d postprocess.d

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

NML_DIR := ../lib
include $(NML_DIR)/nml_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

GEOM_DIR := ../geom
include $(GEOM_DIR)/geom_files.mk

GASDYN_DIR := ../gasdyn
include $(GASDYN_DIR)/gasdyn_files.mk

CEQ_DIR := ../extern/ceq/source
LIBCEQ := $(CEQ_DIR)/libceq.a
include $(CEQ_DIR)/ceq_files.mk

NTYPES_DIR := ../ntypes
NTYPES_FILES := $(NTYPES_DIR)/complex.d

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib

# The install destination.
INSTALL_DIR ?= $(HOME)/gdtkinst

# The build destination sits locally for present
BUILD_DIR := ../../build

REPO_DIR := $(shell cd ../../; pwd)
ifneq ("$(wildcard $(REPO_DIR)/.git)","")
    # Repository exists
    REVISION_STRING := $(shell git rev-parse --short HEAD)
else
    REVISION_STRING := unknown
endif

ifeq ($(DMD), dmd)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -w -g -debug -version=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        DFLAGS := -profile -w -g -O -release -boundscheck=off -version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O -release -boundscheck=off -version=flavour_fast
    endif
    OF := -of
    DVERSION := -version=
    DLINKFLAGS := -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif
ifeq ($(DMD), ldmd2)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -w -g -debug -version=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        DFLAGS := -profile -w -O2 -release -inline -boundscheck=off -version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O2 -release -inline -boundscheck=off -version=flavour_fast
    endif
    OF := -of
    DVERSION := -version=
    DLINKFLAGS := -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif
ifeq ($(DMD), ldc2)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -w -g -d-debug -d-version=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        # -fprofile-generate will result in profraw files being written
        # that may be viewed, showing the top 10 functions with internal block counts
        # llvm-profdata show -text -topn=10 <profraw-file>
        DFLAGS := -fprofile-generate -g -w -O2 -release -enable-inlining -boundscheck=off -d-version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O2 -release -enable-inlining -boundscheck=off -d-version=flavour_fast
        ifeq ($(PLATFORM), linux)
	    # Full link-time optimization does not play nicely on macOS
	    DFLAGS += -flto=full
        endif
    endif
    OF := -of=
    DVERSION := -d-version=
    DLINKFLAGS := -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif
# DIP1008 allows throwing of exceptions in @nogc code. Appeared in 2.079.0, 2018-03-01.
# This rules out the use of gdc for compiling the code.
DFLAGS += -dip1008 -preview=in
DFLAGS += -I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GRID_DIR) -I$(GZIP_DIR)

ifeq ($(WITH_FPE),1)
    DFLAGS += $(DVERSION)enable_fp_exceptions
endif

default: $(PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "L1d4 simulation code built."

install: $(PROGRAMS) prep-gas prep-chem prep-kinetics chemkin2eilmer libgas
	- mkdir -p $(INSTALL_DIR)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share
	cp $(PROGRAMS) $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/bin/* $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/lib/lpeg.so $(BUILD_DIR)/lib
	cp -r ../lib/* $(BUILD_DIR)/lib
	cp $(NML_LUA_MODULES) $(BUILD_DIR)/lib
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

clean:
	- rm *.o
	- rm $(PROGRAMS)
	- rm -r $(BUILD_DIR)/*
	- rm main_with_rev_string.d
	- cd $(LUA_DIR); make clean;
	- cd $(GAS_DIR); make clean; rm libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(CEQ_DIR); make clean

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

$(LIBGASF):
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) PLATFORM=$(PLATFORM) libgasf.a

$(LIBCEQ):
	cd $(CEQ_DIR); make

l1d4: main.d $(L1D_FILES) \
	$(GEOM_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		-e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
		main.d > main_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)l1d4 \
		main_with_rev_string.d \
		$(L1D_FILES) \
		$(GEOM_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBCEQ) $(LIBLUA) \
		$(DLINKFLAGS)

l1d4-prep: l1d4-prep.py
	cp l1d4-prep.py l1d4-prep
	chmod +x l1d4-prep

l1d4-xt-viewer: xt-viewer.py
	cp xt-viewer.py l1d4-xt-viewer
	chmod +x l1d4-xt-viewer

prep-gas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) build-prep-gas

libgas:
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) build-libgas

prep-chem:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-prep-chem

prep-kinetics:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-prep-kinetics

chemkin2eilmer:
	cd $(KINETICS_DIR); make BUILD_DIR=$(BUILD_DIR) build-chemkin2eilmer
