# makefile for the Numerical Types package
#
# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk
NUMBER_FILE := $(filter %number.d,$(NM_FILES))

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NTYPES_FILES := complex.d

ifeq ($(DMD), dmd)
    # DFLAGS := -w
    DFLAGS := -w -O -release -inline -boundscheck=off
    OF := -of
    DLINKFLAGS := -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -O -release -inline -boundscheck=off
    OF := -of
    DLINKFLAGS := -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif

.PHONY: test
test: test_runner
	./test_runner --module ntypes

# We need to pull the test_runner.d to the front to ensure no linker errors
test_runner: $(NUMBER_FILE) $(NTYPES_FILES) $(UTIL_FILES)
	@echo "Compiling test runner..."
	$(DMD) -of$@ -debug -g -w -dip1008 --unittest -version=complex_numbers $(DLINKFLAGS) \
		$(filter %/test_runner.d,$^) $(filter-out %/test_runner.d,$^)

clean:
	- rm test_runner
