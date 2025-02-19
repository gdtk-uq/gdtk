# makefile for slf 
# NNG 2024-05-15

DMD ?= ldc2

# FLAVOUR options are debug, fast, profile
# Flags for each compiler will be determined on this option.
FLAVOUR ?= debug
PLATFORM ?= linux

PROGRAMS := slf

SLF_FILES := io.d linalg.d misc.d

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk
LIBGASF := $(GAS_DIR)/libgasf.a

KINETICS_DIR := ../kinetics
include $(KINETICS_DIR)/kinetics_files.mk

CEQ_DIR := ../extern/ceq/source
LIBCEQ := $(CEQ_DIR)/libceq.a
include $(CEQ_DIR)/ceq_files.mk

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := $(LUA_DIR)/install/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib

DYAML_DIR := $(DGD_REPO)/src/extern/D-YAML/source/dyaml
include $(DYAML_DIR)/dyaml_files.mk

TINYENDIAN_DIR := $(DGD_REPO)/src/extern/tinyendian/source
include $(TINYENDIAN_DIR)/tinyendian_files.mk

# The install destination.
INSTALL_DIR ?= $(HOME)/gdtkinst

# The build destination sits locally for present
BUILD_DIR := ../../build
BUILD_DATE := $(shell date)
REPO_DIR := $(shell cd ../../; pwd)
ifneq ("$(wildcard $(REPO_DIR)/.git)","")
    # Repository exists
    REVISION_STRING := $(shell git rev-parse --short HEAD)
    REVISION_DATE := $(shell git log -1 --format=%cd)
else
    REVISION_STRING := unknown
    REVISION_DATE := unknown
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
    PIC := -fPIC
    OF := -of
    DVERSION := -version=
    DLINKFLAGS :=
    DLINKFLAGS := $(DLINKFLAGS) -L-ldl
endif
ifeq ($(DMD), ldmd2)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -w -g -debug -version=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        DFLAGS := -profile -w -O -release -inline -boundscheck=off -version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O -release -inline -boundscheck=off -version=flavour_fast
    endif
    PIC := -fPIC
    OF := -of
    DVERSION := -version=
    DLINKFLAGS :=
    DLINKFLAGS := $(DLINKFLAGS) -L-ldl
endif
ifeq ($(DMD), ldc2)
    ifeq ($(FLAVOUR), debug)
        DFLAGS := -w -g -link-defaultlib-debug -L-export-dynamic -d-debug -d-version=flavour_debug
    endif
    ifeq ($(FLAVOUR), profile)
        # -fprofile-generate will result in profraw files being written
        # that may be viewed, showing the top 10 functions with internal block counts
        # llvm-profdata show -text -topn=10 <profraw-file>
        DFLAGS := -fprofile-generate -g -w -O -release -enable-inlining -boundscheck=off -d-version=flavour_profile
    endif
    ifeq ($(FLAVOUR), fast)
        DFLAGS := -w -g -O3 -mcpu=native -release -enable-inlining -boundscheck=off -d-version=flavour_fast -flto=full
    endif
    PIC := --relocation-model=pic
    OF := -of=
    DVERSION := -d-version=
    DLINKFLAGS :=
    #ifeq ($(FLAVOUR), profile)
    #    DLINKFLAGS := $(DLINKFLAGS) -Wl,-fprofile-generate
    #endif
    DLINKFLAGS := $(DLINKFLAGS) -L-ldl
endif
# DIP1008 allows throwing of exceptions in @nogc code. Appeared in 2.079.0, 2018-03-01.
# This rules out the use of gdc for compiling the code.
# gdc-8 built on 2.069.2 with updates to 2016-01-15
# gdc-9 built on 2.076.0 which appeared 2017-09-01
# gdc-10 will back-port static foreach, which we also use.
# Expect gdc-10 to be out in May 2020.
DFLAGS += -dip1008 -preview=in
DFLAGS += -I.. -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GRID_DIR) -I$(GZIP_DIR)

ifeq ($(WITH_FPE),1)
    DFLAGS += $(DVERSION)enable_fp_exceptions
endif

default: $(PROGRAMS) $(SLF_FILES)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "slf code built."

install: $(PROGRAMS)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/lib/gdtk
	- cp slf $(BUILD_DIR)/bin
	- cp libslf.so $(BUILD_DIR)/lib
	- cp slf.py $(BUILD_DIR)/lib/gdtk
	cp $(LUA_DIR)/install/bin/* $(BUILD_DIR)/bin
	cp $(LUA_DIR)/install/lib/lpeg.so $(BUILD_DIR)/lib
	cp -r ../lib/* $(BUILD_DIR)/lib
	@echo "Installing to $(INSTALL_DIR)"
	- mkdir -p $(INSTALL_DIR)
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

clean:
	- rm *.o *.so
	- rm $(PROGRAMS)
	- rm -r $(BUILD_DIR)/*
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a
	- cd $(KINETICS_DIR); make clean
	- cd $(CEQ_DIR); make clean
	- cd $(DYAML_DIR); rm -rf *.o
	- cd $(TINYENDIAN_DIR); rm -rf *.o

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local

$(LIBGASF):
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) libgasf.a

$(LIBCEQ):
	cd $(CEQ_DIR); make

slf: slf.d $(SLF_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES)
	$(DMD) $(DFLAGS) $(OF)slf \
		$(DVERSION)complex_numbers \
		slf.d \
		$(SLF_FILES) $(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		$(DLINKFLAGS)

libslf: cwrap.d slf.d $(SLF_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES) \
	$(GAS_FILES) $(CEQ_SRC_FILES) $(LIBCEQ) $(LIBGASF) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(UTIL_FILES)
	$(DMD) -c cwrap.d $(PIC) $(DFLAGS) -op \
		$(DVERSION)complex_numbers \
		slf.d \
		$(SLF_FILES) $(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA)
	$(DMD) $(OF)libslf.so cwrap.o $(PIC) $(DFLAGS) -op -shared \
		$(DVERSION)complex_numbers \
		slf.d \
		$(SLF_FILES) $(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(GAS_FILES) $(CEQ_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(LIBCEQ) $(LIBGASF) $(LIBLUA) \
		-L-L=$(INSTALL_DIR)/lib $(DLINKFLAGS)
